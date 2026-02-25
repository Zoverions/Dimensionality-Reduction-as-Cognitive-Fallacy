"""
Resolution Ladder (Fast) — N=1k and N=5k
=========================================
Optimized for sandbox execution:
  - 10 scales (not 20)
  - 3 bootstraps (not 5)
  - Vectorized pair processing
  - Skips N=15k (infeasible in sandbox; ~days of compute)

Key question: Does the EI_norm behavior change qualitatively between N=1k and N=5k?
"""

import os
import numpy as np
from scipy.spatial import cKDTree
from scipy.stats import entropy
from scipy.sparse import lil_matrix
import matplotlib.pyplot as plt
import time

# --- CONFIGURATION ---
R_VALS      = np.linspace(2.0, 60.0, 10)  # 10 scales
BOX_SIZE    = 205.0 / 0.6774   # TNG300-1 (cMpc/h)
N_BOOTSTRAP = 3
EPS         = 1e-12
DATA_DIR    = '/home/ubuntu/Dimensionality-Reduction-as-Cognitive-Fallacy'


def get_stationary_distribution(W_csr, tol=1e-9, max_iter=200):
    """Exact stationary distribution via Power Iteration: pi W = pi."""
    N = W_csr.shape[0]
    v = np.ones(N) / N
    for _ in range(max_iter):
        v_next = v @ W_csr
        s = v_next.sum()
        if s > 0:
            v_next /= s
        if np.linalg.norm(v_next - v) < tol:
            return v_next
        v = v_next
    return v


def get_exact_sdei(pos, mass, tree, r_link):
    """SDEI with exact stationarity and vectorized graph construction."""
    N = len(pos)
    pairs = tree.query_pairs(r_link, output_type='ndarray')
    
    if len(pairs) == 0:
        return 0.0
    
    i_idxs, j_idxs = pairs[:, 0], pairs[:, 1]
    
    # Periodic boundary conditions
    delta = pos[i_idxs] - pos[j_idxs]
    delta = np.where(delta >  BOX_SIZE/2, delta - BOX_SIZE, delta)
    delta = np.where(delta < -BOX_SIZE/2, delta + BOX_SIZE, delta)
    dists = np.linalg.norm(delta, axis=1)
    
    mask = dists > 0
    i_idxs, j_idxs, dists = i_idxs[mask], j_idxs[mask], dists[mask]
    
    if len(i_idxs) == 0:
        return 0.0
    
    # Gravity weights (asymmetric)
    w_ij = mass[j_idxs] / (dists**2)
    w_ji = mass[i_idxs] / (dists**2)
    
    # Build sparse adjacency
    A = lil_matrix((N, N), dtype=np.float64)
    for k in range(len(i_idxs)):
        A[i_idxs[k], j_idxs[k]] = w_ij[k]
        A[j_idxs[k], i_idxs[k]] = w_ji[k]
    A = A.tocsr()
    
    # Row-stochastic transition matrix
    row_sums = np.array(A.sum(axis=1)).flatten()
    row_sums[row_sums == 0] = 1.0
    W = A.multiply((1.0 / row_sums)[:, np.newaxis]).tocsr()
    
    # EI = H(<W_i>) - <H(W_i)>
    W_mean = np.array(W.mean(axis=0)).flatten()
    H_eff = entropy(W_mean + EPS)
    
    row_entropies = np.zeros(N)
    for i in range(N):
        start, end = W.indptr[i], W.indptr[i+1]
        if end > start:
            row_entropies[i] = entropy(W.data[start:end] + EPS)
    H_noise = np.mean(row_entropies)
    
    EI = H_eff - H_noise
    
    # Normalization: H(pi) from exact stationary distribution
    pi = get_stationary_distribution(W)
    H_rw = entropy(pi + EPS)
    
    return EI / H_rw if H_rw > 0 else 0.0


def compute_ei_curve(pos, mass, label=""):
    """Compute SDEI at each scale."""
    tree = cKDTree(pos, boxsize=BOX_SIZE)
    curve = []
    for i, r in enumerate(R_VALS):
        t0 = time.time()
        val = get_exact_sdei(pos, mass, tree, r)
        dt = time.time() - t0
        curve.append(val)
        print(f"    {label} Scale {i+1}/{len(R_VALS)} ({r:.1f} cMpc/h): EI={val:.4f} [{dt:.1f}s]")
    return np.array(curve)


def run_ladder():
    """Run N=1k and N=5k resolution steps."""
    # Load 15k data (take subsets)
    p0_full = np.load(os.path.join(DATA_DIR, "tng300_z=0_pos_15k.npy"))
    m0_full = np.load(os.path.join(DATA_DIR, "tng300_z=0_mass_15k.npy"))
    p1_full = np.load(os.path.join(DATA_DIR, "tng300_z=1_pos_15k.npy"))
    m1_full = np.load(os.path.join(DATA_DIR, "tng300_z=1_mass_15k.npy"))
    
    print(f"Loaded z=0: {p0_full.shape}, z=1: {p1_full.shape}")
    print(f"Box size: {BOX_SIZE:.2f} cMpc/h")
    
    subsets = [1000, 5000]
    results = {}
    
    for n in subsets:
        n_actual = min(n, len(p0_full), len(p1_full))
        mean_sep = BOX_SIZE / (n_actual ** (1/3))
        
        print(f"\n{'='*60}")
        print(f"  N={n_actual}  |  Mean separation: {mean_sep:.1f} cMpc/h")
        print(f"{'='*60}")
        
        p0, m0 = p0_full[:n_actual], m0_full[:n_actual]
        p1, m1 = p1_full[:n_actual], m1_full[:n_actual]
        
        # Real curves
        print(f"  Real z=0:")
        ei0_real = compute_ei_curve(p0, m0, "z=0")
        print(f"  Real z=1:")
        ei1_real = compute_ei_curve(p1, m1, "z=1")
        delta_real = ei0_real - ei1_real
        
        # Mass-Shuffled Null
        null_deltas = []
        for b in range(N_BOOTSTRAP):
            print(f"  Bootstrap {b+1}/{N_BOOTSTRAP}:")
            m0_shuff = np.random.permutation(m0)
            m1_shuff = np.random.permutation(m1)
            ei0_null = compute_ei_curve(p0, m0_shuff, f"null-z0-b{b+1}")
            ei1_null = compute_ei_curve(p1, m1_shuff, f"null-z1-b{b+1}")
            null_deltas.append(ei0_null - ei1_null)
        
        delta_null_mean = np.mean(null_deltas, axis=0)
        delta_null_std  = np.std(null_deltas, axis=0)
        gain = delta_real - delta_null_mean
        
        results[n_actual] = {
            'gain': gain,
            'delta_real': delta_real,
            'delta_null_mean': delta_null_mean,
            'delta_null_std': delta_null_std,
            'ei0_real': ei0_real,
            'ei1_real': ei1_real,
            'mean_sep': mean_sep,
        }
        
        pk_idx = np.argmax(np.abs(gain))
        sig_mask = np.abs(gain) > 2.0 * (delta_null_std + EPS)
        print(f"\n  N={n_actual} RESULT:")
        print(f"    Peak Gain: {gain[pk_idx]:.6f} at {R_VALS[pk_idx]:.1f} cMpc/h")
        print(f"    Scales > 2σ: {int(sig_mask.sum())}/{len(R_VALS)}")
        print(f"    z=0 EI range: [{ei0_real.min():.4f}, {ei0_real.max():.4f}]")
        print(f"    z=1 EI range: [{ei1_real.min():.4f}, {ei1_real.max():.4f}]")
    
    return results


def plot_results(results, out_path):
    """3-panel figure: Gain, EI_norm z=0, EI_norm z=1."""
    fig, axes = plt.subplots(1, 3, figsize=(20, 6))
    colors = {1000: '#1f77b4', 5000: '#ff7f0e', 15000: '#2ca02c'}
    
    # Panel 1: Topological Gain
    ax = axes[0]
    for n, res in results.items():
        c = colors.get(n, 'gray')
        ax.plot(R_VALS, res['gain'], '-o', color=c, linewidth=2, markersize=4,
                label=f'N={n} (sep={res["mean_sep"]:.0f} cMpc/h)')
        ax.fill_between(R_VALS,
                        -2 * res['delta_null_std'],
                        +2 * res['delta_null_std'],
                        color=c, alpha=0.1)
    ax.axhline(0, color='black', linestyle='--', alpha=0.8, label='Null')
    ax.set_xlabel('Linking Length (cMpc/h)', fontsize=11)
    ax.set_ylabel('Topological Gain', fontsize=11)
    ax.set_title('Double Difference: Does the Signal Emerge?', fontsize=12)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    
    # Panel 2: Absolute EI_norm at z=0
    ax2 = axes[1]
    for n, res in results.items():
        c = colors.get(n, 'gray')
        ax2.plot(R_VALS, res['ei0_real'], '-o', color=c, linewidth=2, markersize=4,
                 label=f'N={n}')
    ax2.set_xlabel('Linking Length (cMpc/h)', fontsize=11)
    ax2.set_ylabel('EI_norm', fontsize=11)
    ax2.set_title('Absolute EI_norm at z=0 (Mature Web)', fontsize=12)
    ax2.legend(fontsize=9)
    ax2.grid(True, alpha=0.3)
    
    # Panel 3: Absolute EI_norm at z=1
    ax3 = axes[2]
    for n, res in results.items():
        c = colors.get(n, 'gray')
        ax3.plot(R_VALS, res['ei1_real'], '-o', color=c, linewidth=2, markersize=4,
                 label=f'N={n}')
    ax3.set_xlabel('Linking Length (cMpc/h)', fontsize=11)
    ax3.set_ylabel('EI_norm', fontsize=11)
    ax3.set_title('Absolute EI_norm at z=1 (Early Web)', fontsize=12)
    ax3.legend(fontsize=9)
    ax3.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    print(f"\nFigure saved to {out_path}")
    
    # Summary table
    print("\n" + "="*70)
    print(f"{'N':>8s}  {'MeanSep':>8s}  {'PeakGain':>10s}  {'PeakScale':>10s}  {'Sig>2σ':>8s}")
    print("="*70)
    for n, res in results.items():
        gain = res['gain']
        pk = gain[np.argmax(np.abs(gain))]
        ps = R_VALS[np.argmax(np.abs(gain))]
        ns = int(np.sum(np.abs(gain) > 2 * (res['delta_null_std'] + EPS)))
        print(f"{n:>8d}  {res['mean_sep']:>8.1f}  {pk:>10.6f}  {ps:>10.1f}  {ns:>4d}/{len(R_VALS)}")
    print("="*70)


if __name__ == "__main__":
    np.random.seed(42)
    results = run_ladder()
    out_path = os.path.join(DATA_DIR, 'resolution_ladder.png')
    plot_results(results, out_path)
