"""
Spatial Shuffle Diagnostic
==========================
Tests whether the Double-Difference null result is a resolution artifact
or a genuine property of the operator.

Three conditions compared at each redshift:
  1. Real:           Real positions + Real masses
  2. Mass-Shuffled:  Real positions + Shuffled masses  (existing null)
  3. Spatial-Shuffled: Shuffled positions + Real masses (NEW diagnostic)

The spatial shuffle randomizes halo positions uniformly within the box
while preserving the exact mass function. This destroys all spatial
clustering while keeping the mass distribution identical.

Interpretation:
  - If Spatial-Shuffled EI_norm ≈ Real EI_norm → the operator is
    insensitive to spatial structure at this N (genuine null)
  - If Spatial-Shuffled EI_norm ≠ Real EI_norm → the operator DOES
    detect spatial structure, and the mass-shuffled null result means
    the mass-position correlation is not yet detectable (resolution floor)
"""

import os
import numpy as np
from scipy.spatial import cKDTree
from scipy.stats import entropy
from scipy.sparse import lil_matrix
import matplotlib.pyplot as plt

# --- CONFIGURATION ---
BOX_SIZE     = 205.0 / 0.6774       # TNG300-1 box in cMpc/h
N_SCALES     = 15
R_VALS       = np.linspace(2.0, 60.0, N_SCALES)
EPS          = 1e-12
N_SPATIAL_BOOT = 5   # Number of spatial shuffle realizations
DATA_DIR     = '/home/ubuntu/Dimensionality-Reduction-as-Cognitive-Fallacy'

# --- EI OPERATOR (same as double_difference.py) ---
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

    A[i_idxs, j_idxs] = mass[j_idxs] / (dists**2)
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


def compute_ei_curve(pos, mass):
    tree = cKDTree(pos, boxsize=BOX_SIZE)
    return np.array([get_normalized_ei(pos, mass, tree, r) for r in R_VALS])


def run_diagnostic(label, pos, mass):
    N = len(pos)
    print(f"\n=== {label} (N={N}) ===")

    # 1. Real
    print(f"  Computing Real EI_norm...")
    ei_real = compute_ei_curve(pos, mass)

    # 2. Mass-Shuffled (existing null)
    print(f"  Computing Mass-Shuffled null...")
    ei_mass_shuffled = compute_ei_curve(pos, np.random.permutation(mass))

    # 3. Spatial-Shuffled (NEW diagnostic)
    print(f"  Computing {N_SPATIAL_BOOT} Spatial-Shuffled realizations...")
    ei_spatial_shuffles = []
    for b in range(N_SPATIAL_BOOT):
        pos_rand = np.random.rand(N, 3) * BOX_SIZE
        ei_spatial_shuffles.append(compute_ei_curve(pos_rand, mass))
        print(f"    Spatial shuffle {b+1}/{N_SPATIAL_BOOT} done")

    ei_spatial_mean = np.mean(ei_spatial_shuffles, axis=0)
    ei_spatial_std  = np.std(ei_spatial_shuffles, axis=0)

    return {
        'ei_real':          ei_real,
        'ei_mass_shuffled': ei_mass_shuffled,
        'ei_spatial_mean':  ei_spatial_mean,
        'ei_spatial_std':   ei_spatial_std,
        'N':                N,
    }


def plot_diagnostic(res_z0, res_z1, out_path):
    fig, axes = plt.subplots(1, 2, figsize=(16, 6), sharey=True)

    for ax, res, label in [(axes[0], res_z0, 'z=0 (Mature Web)'),
                           (axes[1], res_z1, 'z=1 (Early Web)')]:
        ax.plot(R_VALS, res['ei_real'], 'b-o', linewidth=2, markersize=4,
                label='Real (pos + mass)')
        ax.plot(R_VALS, res['ei_mass_shuffled'], 'g--s', linewidth=1.5,
                markersize=4, alpha=0.8, label='Mass-Shuffled (real pos)')
        ax.plot(R_VALS, res['ei_spatial_mean'], 'r-^', linewidth=2,
                markersize=4, label='Spatial-Shuffled (random pos)')
        ax.fill_between(R_VALS,
                        res['ei_spatial_mean'] - 2 * res['ei_spatial_std'],
                        res['ei_spatial_mean'] + 2 * res['ei_spatial_std'],
                        color='red', alpha=0.1, label='Spatial-Shuffled +/-2s')

        ax.set_xlabel('Linking Length (cMpc/h)', fontsize=12)
        ax.set_title(f'{label} (N={res["N"]})', fontsize=12)
        ax.legend(fontsize=8, loc='best')
        ax.grid(True, alpha=0.3)

    axes[0].set_ylabel('EI_norm', fontsize=12)
    fig.suptitle('Spatial Shuffle Diagnostic: Does the Operator Detect Structure?',
                 fontsize=13, fontweight='bold')
    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    print(f"\nFigure saved to {out_path}")

    # Print interpretation
    for label, res in [('z=0', res_z0), ('z=1', res_z1)]:
        real_peak = res['ei_real'].max()
        spatial_peak = res['ei_spatial_mean'].max()
        ratio = real_peak / spatial_peak if spatial_peak > 0 else float('inf')
        print(f"\n  {label}:")
        print(f"    Real peak EI_norm:           {real_peak:.4f}")
        print(f"    Spatial-Shuffled peak:        {spatial_peak:.4f}")
        print(f"    Ratio (Real / Spatial):       {ratio:.4f}")
        if ratio < 0.5:
            print(f"    --> OPERATOR DETECTS STRUCTURE (Real << Spatial)")
        elif ratio > 1.5:
            print(f"    --> ANOMALY: Real > Spatial (unexpected)")
        else:
            print(f"    --> INCONCLUSIVE: Real ~ Spatial")


if __name__ == "__main__":
    np.random.seed(42)

    pos_z0  = np.load(os.path.join(DATA_DIR, "tng300_z=0_pos.npy"))
    mass_z0 = np.load(os.path.join(DATA_DIR, "tng300_z=0_mass.npy"))
    pos_z1  = np.load(os.path.join(DATA_DIR, "tng300_z=1_pos.npy"))
    mass_z1 = np.load(os.path.join(DATA_DIR, "tng300_z=1_mass.npy"))

    res_z0 = run_diagnostic("z=0", pos_z0, mass_z0)
    res_z1 = run_diagnostic("z=1", pos_z1, mass_z1)

    plot_diagnostic(res_z0, res_z1,
                    os.path.join(DATA_DIR, 'spatial_shuffle_diagnostic.png'))
