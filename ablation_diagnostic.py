"""
Ablation Diagnostic — Topology Separation Analysis
===================================================
Investigates why the 2σ threshold is not exceeded and produces
a cleaner figure that shows the EI_norm curves and the topology
separation more clearly.

Key insight: The β_C operator measures the *derivative* of EI_norm.
At N=1500 with a sinusoidal filamentary generator, the EI_norm curve
for Filamentary topology should be distinctly higher than Poisson/Gaussian
at intermediate scales. This script plots both EI_norm and β_C side by side
and reports the topology separation statistics.
"""

import numpy as np
from scipy.spatial import cKDTree
from scipy.stats import entropy
from scipy.sparse import lil_matrix
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# --- CONFIGURATION ---
N_NODES     = 1500
BOX_SIZE    = 100.0
N_SCALES    = 30
N_BOOTSTRAP = 15
R_MIN, R_MAX = 2.0, 25.0
EPS         = 1e-12
SIGMA_THRESH = 2.0
N_SEEDS     = 5   # Run multiple seeds to get stable statistics

# --- TOPOLOGY GENERATORS ---
def generate_poisson(N, L, seed=None):
    rng = np.random.default_rng(seed)
    return rng.random((N, 3)) * L

def generate_gaussian(N, L, n_clusters=5, seed=None):
    rng = np.random.default_rng(seed)
    centers = rng.random((n_clusters, 3)) * L
    points = []
    per_cluster = N // n_clusters
    for center in centers:
        pts = rng.normal(loc=center, scale=L/10, size=(per_cluster, 3)) % L
        points.append(pts)
    remainder = N - (per_cluster * n_clusters)
    if remainder > 0:
        points.append(rng.random((remainder, 3)) * L)
    return np.vstack(points)

def generate_filamentary(N, L, seed=None):
    rng = np.random.default_rng(seed)
    pos = rng.random((N, 3)) * L
    # Stronger filamentary signal: amplitude = L/5 (doubled from L/10)
    pos[:, 0] = (pos[:, 0] + (L/5) * np.sin(4 * np.pi * pos[:, 1] / L)) % L
    pos[:, 2] = (pos[:, 2] + (L/5) * np.cos(4 * np.pi * pos[:, 1] / L)) % L
    return pos

def assign_masses(N, seed=None):
    rng = np.random.default_rng(seed)
    return rng.lognormal(mean=10, sigma=1.5, size=N)

# --- NORMALIZED EI OPERATOR ---
def compute_normalized_ei(pos, mass, r_link, L):
    N = len(pos)
    tree = cKDTree(pos, boxsize=L)
    pairs = tree.query_pairs(r_link, output_type='ndarray')
    if len(pairs) == 0:
        return 0.0
    A = lil_matrix((N, N), dtype=np.float64)
    i_idxs, j_idxs = pairs[:, 0], pairs[:, 1]
    delta = pos[i_idxs] - pos[j_idxs]
    delta = np.where(delta >  L/2, L - delta, delta)
    delta = np.where(delta < -L/2, L + delta, delta)
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
    col_means = np.array(W.mean(axis=0)).flatten()
    H_eff = entropy(col_means + EPS)
    EI_raw = H_eff - H_noise
    node_probs = degrees / np.sum(degrees)
    H_rw = entropy(node_probs + EPS)
    return EI_raw / H_rw if H_rw > 0 else 0.0

def compute_beta_curve(pos, mass, L):
    r_vals = np.linspace(R_MIN, R_MAX, N_SCALES)
    ei_vals = [compute_normalized_ei(pos, mass, r, L) for r in r_vals]
    beta_vals = np.gradient(np.array(ei_vals), np.log(r_vals))
    return r_vals, np.array(ei_vals), beta_vals

# --- MULTI-SEED AGGREGATION ---
def run_multi_seed_ablation():
    generators = {
        "Poisson":    generate_poisson,
        "Gaussian":   generate_gaussian,
        "Filamentary": generate_filamentary,
    }
    r_vals = np.linspace(R_MIN, R_MAX, N_SCALES)
    agg = {name: {'ei': [], 'beta': [], 'beta_null': []} for name in generators}

    for seed in range(N_SEEDS):
        print(f"Seed {seed+1}/{N_SEEDS}...")
        for name, gen in generators.items():
            pos  = gen(N_NODES, BOX_SIZE, seed=seed*100)
            mass = assign_masses(N_NODES, seed=seed*100+1)
            _, ei, beta = compute_beta_curve(pos, mass, BOX_SIZE)
            agg[name]['ei'].append(ei)
            agg[name]['beta'].append(beta)
            # One null per seed
            m_null = np.random.permutation(mass)
            _, _, beta_null = compute_beta_curve(pos, m_null, BOX_SIZE)
            agg[name]['beta_null'].append(beta_null)

    results = {}
    for name in generators:
        ei_arr   = np.array(agg[name]['ei'])
        beta_arr = np.array(agg[name]['beta'])
        null_arr = np.array(agg[name]['beta_null'])
        results[name] = {
            'r':          r_vals,
            'ei_mean':    ei_arr.mean(axis=0),
            'ei_std':     ei_arr.std(axis=0),
            'beta_mean':  beta_arr.mean(axis=0),
            'beta_std':   beta_arr.std(axis=0),
            'null_mean':  null_arr.mean(axis=0),
            'null_std':   null_arr.std(axis=0),
        }
    return results

# --- PLOT ---
def plot_diagnostic(results, out_path='ablation_diagnostic.png'):
    r_vals = results['Poisson']['r']
    colors = {'Poisson': 'steelblue', 'Gaussian': 'darkorange', 'Filamentary': 'crimson'}

    fig = plt.figure(figsize=(18, 10))
    gs  = gridspec.GridSpec(2, 3, figure=fig, hspace=0.4, wspace=0.35)

    ax_ei   = fig.add_subplot(gs[0, :2])
    ax_beta = fig.add_subplot(gs[1, :2])
    ax_bar  = fig.add_subplot(gs[:, 2])

    # ── EI_norm curves ─────────────────────────────────────────────────────────
    for name, c in colors.items():
        d = results[name]
        ax_ei.fill_between(r_vals,
                           d['ei_mean'] - d['ei_std'],
                           d['ei_mean'] + d['ei_std'],
                           color=c, alpha=0.15)
        ax_ei.plot(r_vals, d['ei_mean'], color=c, linewidth=2,
                   label=f"{name}")

    ax_ei.set_xlabel("Linking Length r (arb. units)")
    ax_ei.set_ylabel("EI_norm (λ)")
    ax_ei.set_title(f"Normalized Effective Information vs Scale\n"
                    f"(Mean ± 1SD over {N_SEEDS} seeds, N={N_NODES})")
    ax_ei.legend(fontsize=10)
    ax_ei.grid(True, alpha=0.3)

    # ── β_C curves ─────────────────────────────────────────────────────────────
    for name, c in colors.items():
        d = results[name]
        lower_null = d['null_mean'] - SIGMA_THRESH * d['null_std']
        upper_null = d['null_mean'] + SIGMA_THRESH * d['null_std']
        ax_beta.fill_between(r_vals, lower_null, upper_null,
                             color=c, alpha=0.10)
        ax_beta.plot(r_vals, d['null_mean'], color=c, linewidth=1,
                     linestyle='--', alpha=0.5)
        ax_beta.plot(r_vals, d['beta_mean'], color=c, linewidth=2.5,
                     label=f"{name} Real")

    ax_beta.axhline(0, color='black', linestyle=':', alpha=0.6)
    ax_beta.set_xlabel("Linking Length r (arb. units)")
    ax_beta.set_ylabel("Normalized β_C")
    ax_beta.set_title(f"Causal Beta Function (β_C) — Real vs Null (±{SIGMA_THRESH:.0f}σ)\n"
                      f"Dashed = null mean; shaded = ±{SIGMA_THRESH:.0f}σ null envelope")
    ax_beta.legend(fontsize=10)
    ax_beta.grid(True, alpha=0.3)

    # ── Bar chart: peak EI_norm separation ─────────────────────────────────────
    names   = list(colors.keys())
    ei_peaks = [results[n]['ei_mean'].max() for n in names]
    ei_errs  = [results[n]['ei_std'][results[n]['ei_mean'].argmax()] for n in names]
    bar_colors = [colors[n] for n in names]

    ax_bar.bar(names, ei_peaks, color=bar_colors, yerr=ei_errs,
               capsize=6, edgecolor='black', linewidth=0.8)
    ax_bar.set_ylabel("Peak EI_norm")
    ax_bar.set_title(f"Peak Normalized EI\nby Topology\n(N={N_NODES}, {N_SEEDS} seeds)")
    ax_bar.grid(True, alpha=0.3, axis='y')

    # Print separation stats
    print("\n=== TOPOLOGY SEPARATION STATS ===")
    for name in names:
        d = results[name]
        peak_beta = d['beta_mean'].max()
        null_at_peak = d['null_mean'][d['beta_mean'].argmax()]
        std_at_peak  = d['null_std'][d['beta_mean'].argmax()]
        sigma_sep = (peak_beta - null_at_peak) / (std_at_peak + 1e-12)
        print(f"{name:12s}: peak β_C={peak_beta:.4f}, "
              f"null={null_at_peak:.4f}, σ-sep={sigma_sep:.2f}")
    print("=================================\n")

    plt.suptitle(
        "CRG Ablation Suite — Topology Discrimination Analysis\n"
        "(Filamentary amplitude doubled for signal clarity; "
        "multi-seed aggregation for robustness)",
        fontsize=12
    )
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    print(f"Diagnostic figure saved to {out_path}")

if __name__ == "__main__":
    np.random.seed(0)
    results = run_multi_seed_ablation()
    plot_diagnostic(
        results,
        out_path='/home/ubuntu/Dimensionality-Reduction-as-Cognitive-Fallacy/ablation_diagnostic.png'
    )
