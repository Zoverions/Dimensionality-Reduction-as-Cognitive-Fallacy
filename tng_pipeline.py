# @title Phase 2: Normalized TNG Temporal CRG Sweep (Matched z=1 vs z=0)
import requests
import numpy as np
from scipy.spatial import cKDTree
from scipy.stats import entropy
from scipy.sparse import lil_matrix, csr_matrix
import matplotlib.pyplot as plt

# --- 1. CONFIGURATION ---
API_KEY = "INSERT_YOUR_API_KEY_HERE"
HEADERS = {"api-key": API_KEY}
BASE_URL = "http://www.tng-project.org/api/TNG100-1/snapshots/"
BOX_SIZE = 75.0 / 0.704
LIMIT = 5000
N_BOOTSTRAP = 20
SNAPSHOTS = {'z=1 (Early Web)': 85, 'z=0 (Mature Web)': 99}
R_VALS = np.linspace(2.0, 20.0, 25)
EPS = 1e-12

def fetch_snapshot_data(snap_num):
    print(f"Fetching {LIMIT} matched halos for Snapshot {snap_num}...")
    pos, mass = [], []
    url = f"{BASE_URL}{snap_num}/subhalos/"
    params = {'limit': min(LIMIT, 1000), 'order_by': '-mass', 'fields': 'pos,mass'}

    while url and len(pos) < LIMIT:
        r = requests.get(url, params=params, headers=HEADERS)
        r.raise_for_status()
        data = r.json()
        for s in data['results']:
            pos.append(s['pos'])
            mass.append(s['mass'])
            if len(pos) >= LIMIT: break
        url = data.get('next')
        params = {}
    return np.array(pos) / 1000.0, np.array(mass)

def get_normalized_effective_information(pos, mass, tree, r_link):
    N = len(pos)
    pairs = tree.query_pairs(r_link, output_type='ndarray')
    if len(pairs) == 0: return 0.0, N

    A = lil_matrix((N, N), dtype=np.float64)
    i_idxs, j_idxs = pairs[:, 0], pairs[:, 1]

    delta = pos[i_idxs] - pos[j_idxs]
    delta = np.where(delta > BOX_SIZE/2, BOX_SIZE - delta, delta)
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

    H_noise = np.mean(row_entropies)
    W_eff = np.array(W.mean(axis=0)).flatten()
    H_eff = entropy(W_eff + EPS)

    EI_raw = H_eff - H_noise

    node_probabilities = degrees / np.sum(degrees)
    H_rw = entropy(node_probabilities + EPS)

    EI_norm = EI_raw / H_rw if H_rw > 0 else 0.0
    effective_support = np.sum(degrees > 1.0) # Track cardinality

    return EI_norm, effective_support

def run_temporal_pipeline():
    if API_KEY == "INSERT_YOUR_API_KEY_HERE":
        print("Error: Please insert a valid IllustrisTNG API Key in the script.")
        return None

    results = {}
    for label, snap in SNAPSHOTS.items():
        print(f"\n--- Processing {label} ---")
        pos_real, mass_real = fetch_snapshot_data(snap)
        tree = cKDTree(pos_real, boxsize=BOX_SIZE)

        ei_real = []
        for r in R_VALS:
            ei_n, support = get_normalized_effective_information(pos_real, mass_real, tree, r)
            ei_real.append(ei_n)

        beta_real = np.gradient(ei_real, np.log(R_VALS))

        beta_rands = []
        for b in range(N_BOOTSTRAP):
            m_rand = np.copy(mass_real); np.random.shuffle(m_rand)
            ei_r = [get_normalized_effective_information(pos_real, m_rand, tree, r)[0] for r in R_VALS]
            beta_rands.append(np.gradient(ei_r, np.log(R_VALS)))

        results[label] = {
            'beta_real': beta_real,
            'mu_null': np.mean(beta_rands, axis=0),
            'std_null': np.std(beta_rands, axis=0)
        }
    return results

def plot_temporal_evolution(results):
    if results is None: return

    fig, ax = plt.subplots(figsize=(10, 6))
    colors = {'z=1 (Early Web)': 'orange', 'z=0 (Mature Web)': 'blue'}

    for label, data in results.items():
        c = colors[label]
        ax.fill_between(R_VALS, data['mu_null'] - data['std_null'], data['mu_null'] + data['std_null'], color=c, alpha=0.15, label=f"{label} Control (±1 SD)")
        ax.plot(R_VALS, data['beta_real'], color=c, linewidth=2.5, marker='o', markersize=4, label=f"{label} Real Topology")

    ax.axvspan(7.0, 13.0, color='gray', alpha=0.1, label='Typical Sponge Regime')
    ax.axhline(0, color='black', linestyle=':', alpha=0.8)
    ax.set_xlabel('Spatial Scale λ (Mpc)', fontsize=12)
    ax.set_ylabel('Normalized Causal Beta (β_C)', fontsize=12)
    ax.set_title('Evolution of Normalized Causal Resurgence over Cosmic Time', fontsize=14, pad=15)
    ax.legend(loc='upper right', fontsize=10)
    ax.grid(True, alpha=0.3)

    plt.savefig('temporal_evolution.png')
    print("Temporal evolution plot saved to temporal_evolution.png")

if __name__ == "__main__":
    results = run_temporal_pipeline()
    plot_temporal_evolution(results)
