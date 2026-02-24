"""
TNG Temporal CRG Pipeline — Optimized
======================================
Fetches real dark matter subhalo data from IllustrisTNG TNG100-1 via the
public REST API and runs the normalized CRG sweep across two cosmic epochs.

Optimized for sandbox execution:
  - N=500 most massive subhalos per snapshot
  - 10 bootstrap nulls
  - 15 scale points
  - Caches fetched data to .npz files for re-runs
"""

import os
import requests
import numpy as np
from scipy.spatial import cKDTree
from scipy.stats import entropy
from scipy.sparse import lil_matrix
import matplotlib.pyplot as plt
from concurrent.futures import ThreadPoolExecutor, as_completed
import time

# --- 1. CONFIGURATION ---
API_KEY      = "9d963b728d2ffbc4778b7b3f5188e742"
HEADERS      = {"api-key": API_KEY}
BASE_URL     = "https://www.tng-project.org/api/TNG100-1/snapshots/"
H_PARAM      = 0.6774
BOX_SIZE     = 75.0 / H_PARAM          # ~110.7 Mpc
LIMIT        = 500                      # Top 500 most massive subhalos
N_BOOTSTRAP  = 10
N_SCALES     = 15
R_VALS       = np.linspace(2.0, 20.0, N_SCALES)
EPS          = 1e-12
SIGMA_THRESH = 2.0
MAX_WORKERS  = 10
SNAPSHOTS    = {'z~0.5 (Early Web)': 85, 'z=0 (Mature Web)': 99}
CACHE_DIR    = '/home/ubuntu/Dimensionality-Reduction-as-Cognitive-Fallacy'

# --- 2. DATA FETCHING WITH CACHING ---
def fetch_top_ids(snap_num, limit):
    ids = []
    url = f"{BASE_URL}{snap_num}/subhalos/"
    params = {'limit': min(limit, 500), 'order_by': '-mass_log_msun'}
    while url and len(ids) < limit:
        r = requests.get(url, params=params, headers=HEADERS, timeout=60)
        r.raise_for_status()
        data = r.json()
        for s in data['results']:
            ids.append(s['id'])
            if len(ids) >= limit:
                break
        url = data.get('next')
        params = {}
    return ids[:limit]

def fetch_subhalo_detail(snap_num, subhalo_id):
    url = f"{BASE_URL}{snap_num}/subhalos/{subhalo_id}/"
    for attempt in range(3):
        try:
            r = requests.get(url, headers=HEADERS, timeout=30)
            r.raise_for_status()
            d = r.json()
            pos = np.array([d['pos_x'], d['pos_y'], d['pos_z']]) / (1000.0 * H_PARAM)
            mass = d['mass']
            return pos, mass
        except Exception:
            time.sleep(1.5 * (attempt + 1))
    return None, None

def fetch_or_load(snap_num):
    cache_file = os.path.join(CACHE_DIR, f'tng_snap{snap_num}_N{LIMIT}.npz')
    if os.path.exists(cache_file):
        print(f"  Loading cached data from {cache_file}")
        data = np.load(cache_file)
        return data['pos'], data['mass']

    print(f"  Fetching top {LIMIT} subhalo IDs for snapshot {snap_num}...")
    ids = fetch_top_ids(snap_num, LIMIT)
    print(f"  Got {len(ids)} IDs. Fetching positions ({MAX_WORKERS} workers)...")

    pos_list, mass_list = [None] * len(ids), [None] * len(ids)
    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as ex:
        future_to_idx = {ex.submit(fetch_subhalo_detail, snap_num, sid): i
                         for i, sid in enumerate(ids)}
        done = 0
        for fut in as_completed(future_to_idx):
            i = future_to_idx[fut]
            pos, mass = fut.result()
            if pos is not None:
                pos_list[i]  = pos
                mass_list[i] = mass
            done += 1
            if done % 100 == 0:
                print(f"    {done}/{len(ids)} fetched...")

    valid = [(p, m) for p, m in zip(pos_list, mass_list) if p is not None]
    print(f"  Successfully fetched {len(valid)}/{len(ids)} subhalos.")
    pos_arr  = np.array([v[0] for v in valid])
    mass_arr = np.array([v[1] for v in valid])

    np.savez(cache_file, pos=pos_arr, mass=mass_arr)
    print(f"  Cached to {cache_file}")
    return pos_arr, mass_arr

# --- 3. NORMALIZED EI OPERATOR ---
def get_normalized_ei(pos, mass, tree, r_link):
    N = len(pos)
    pairs = tree.query_pairs(r_link, output_type='ndarray')
    if len(pairs) == 0:
        return 0.0

    A = lil_matrix((N, N), dtype=np.float64)
    i_idxs, j_idxs = pairs[:, 0], pairs[:, 1]

    delta = pos[i_idxs] - pos[j_idxs]
    delta = np.where(delta >  BOX_SIZE/2, BOX_SIZE - delta, delta)
    delta = np.where(delta < -BOX_SIZE/2, BOX_SIZE + delta, delta)
    dists = np.linalg.norm(delta, axis=1)

    mask = dists > 0
    i_idxs, j_idxs, dists = i_idxs[mask], j_idxs[mask], dists[mask]

    weights = mass[j_idxs] / (dists**2)
    A[i_idxs, j_idxs] = weights
    A[j_idxs, i_idxs] = mass[i_idxs] / (dists**2)

    A = A.tocsr()
    degrees = np.array(A.sum(axis=1)).flatten()
    degrees[degrees == 0] = 1.0
    W = A.multiply((1.0 / degrees)[:, np.newaxis]).tocsr()

    row_entropies = np.zeros(N)
    for i in range(N):
        start, end = W.indptr[i], W.indptr[i+1]
        if end > start:
            row_entropies[i] = entropy(W.data[start:end] + EPS)

    H_noise    = np.mean(row_entropies)
    W_eff      = np.array(W.mean(axis=0)).flatten()
    H_eff      = entropy(W_eff + EPS)
    EI_raw     = H_eff - H_noise
    node_probs = degrees / np.sum(degrees)
    H_rw       = entropy(node_probs + EPS)
    return EI_raw / H_rw if H_rw > 0 else 0.0

# --- 4. BETA CURVE ---
def compute_beta_curve(pos, mass):
    tree = cKDTree(pos, boxsize=BOX_SIZE)
    ei_vals = [get_normalized_ei(pos, mass, tree, r) for r in R_VALS]
    beta    = np.gradient(np.array(ei_vals), np.log(R_VALS))
    return np.array(ei_vals), beta

# --- 5. MAIN PIPELINE ---
def run_temporal_pipeline():
    results = {}
    for label, snap in SNAPSHOTS.items():
        print(f"\n=== Processing {label} (snap {snap}) ===")
        pos_real, mass_real = fetch_or_load(snap)

        print(f"  Computing real beta curve ({N_SCALES} scales, N={len(pos_real)})...")
        ei_real, beta_real = compute_beta_curve(pos_real, mass_real)

        print(f"  Running {N_BOOTSTRAP} bootstrap nulls...")
        beta_nulls = []
        for b in range(N_BOOTSTRAP):
            m_rand = np.random.permutation(mass_real)
            _, beta_null = compute_beta_curve(pos_real, m_rand)
            beta_nulls.append(beta_null)
            if (b + 1) % 5 == 0:
                print(f"    Bootstrap {b+1}/{N_BOOTSTRAP} done")

        results[label] = {
            'ei_real':   ei_real,
            'beta_real': beta_real,
            'mu_null':   np.mean(beta_nulls, axis=0),
            'std_null':  np.std(beta_nulls, axis=0),
            'N':         len(mass_real),
        }
    return results

# --- 6. SIGNIFICANCE ---
def check_significance(beta_real, mu_null, std_null):
    upper    = mu_null + SIGMA_THRESH * std_null
    peak_idx = np.argmax(beta_real)
    peak_sig = beta_real[peak_idx] > upper[peak_idx]
    n_sig    = int(np.sum(beta_real > upper))
    return peak_sig, n_sig, peak_idx

# --- 7. PLOT ---
def plot_temporal_evolution(results, out_path='temporal_evolution.png'):
    if not results:
        return

    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    colors = {'z~0.5 (Early Web)': 'darkorange', 'z=0 (Mature Web)': 'steelblue'}

    ax = axes[0]
    for label, data in results.items():
        c     = colors[label]
        lower = data['mu_null'] - SIGMA_THRESH * data['std_null']
        upper = data['mu_null'] + SIGMA_THRESH * data['std_null']
        peak_sig, n_sig, _ = check_significance(
            data['beta_real'], data['mu_null'], data['std_null'])
        sig_str = f"Peak>{SIGMA_THRESH:.0f}s: {'YES' if peak_sig else 'NO'}"

        ax.fill_between(R_VALS, lower, upper, color=c, alpha=0.15,
                        label=f"{label} Null (+/-{SIGMA_THRESH:.0f}s)")
        ax.plot(R_VALS, data['mu_null'], color=c, linewidth=1,
                linestyle='--', alpha=0.5)
        ax.plot(R_VALS, data['beta_real'], color=c, linewidth=2.5,
                marker='o', markersize=4,
                label=f"{label} Real (N={data['N']}, {sig_str})")

    ax.axvspan(7.0, 13.0, color='gray', alpha=0.08, label='Sponge Topology Window')
    ax.axhline(0, color='black', linestyle=':', alpha=0.8)
    ax.set_xlabel('Spatial Scale (Mpc)', fontsize=12)
    ax.set_ylabel('Normalized Causal Beta (Beta_C)', fontsize=12)
    ax.set_title(f'Temporal Evolution of Normalized Causal Resurgence\n'
                 f'TNG100-1 Real Subhalo Data (+/-{SIGMA_THRESH:.0f}s null)',
                 fontsize=12)
    ax.legend(loc='best', fontsize=8)
    ax.grid(True, alpha=0.3)

    ax2 = axes[1]
    for label, data in results.items():
        c = colors[label]
        ax2.plot(R_VALS, data['ei_real'], color=c, linewidth=2.5,
                 marker='s', markersize=4, label=f"{label}")
    ax2.set_xlabel('Spatial Scale (Mpc)', fontsize=12)
    ax2.set_ylabel('EI_norm', fontsize=12)
    ax2.set_title('Normalized Effective Information vs Scale\n'
                  '(TNG100-1 Real Subhalo Data)', fontsize=12)
    ax2.legend(fontsize=9)
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    print(f"\nFigure saved to {out_path}")

    print("\n=== TEMPORAL PIPELINE RESULTS ===")
    for label, data in results.items():
        peak_sig, n_sig, peak_idx = check_significance(
            data['beta_real'], data['mu_null'], data['std_null'])
        print(f"{label}:")
        print(f"  N subhalos:         {data['N']}")
        print(f"  Peak Beta_C:        {data['beta_real'].max():.4f}")
        print(f"  Peak scale:         {R_VALS[data['beta_real'].argmax()]:.2f} Mpc")
        print(f"  Null mean at peak:  {data['mu_null'][data['beta_real'].argmax()]:.4f}")
        print(f"  Null std at peak:   {data['std_null'][data['beta_real'].argmax()]:.4f}")
        print(f"  Exceeds {SIGMA_THRESH:.0f}s:         {'YES' if peak_sig else 'NO'}")
        print(f"  Scales above {SIGMA_THRESH:.0f}s:    {n_sig}/{N_SCALES}")
    print("==================================")

if __name__ == "__main__":
    np.random.seed(42)
    results = run_temporal_pipeline()
    plot_temporal_evolution(
        results,
        out_path=os.path.join(CACHE_DIR, 'temporal_evolution.png')
    )
