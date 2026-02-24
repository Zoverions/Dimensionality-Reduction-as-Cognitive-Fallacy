"""
Phase 3: Double-Difference Analysis — Secular Evolution of Causal Structure
============================================================================
Tests whether EI_norm evolves from z=1 to z=0 in a way that cannot be
explained by simple clustering growth.

Two nulls:
  1. Mass-Shuffled Null: preserves exact positions, shuffles masses.
     Controls for geometry — isolates mass-topology correlation.
  2. RGG Null: Poisson-distributed points with same density.
     Controls for clustering — isolates spatial structure.

The Double Difference:
  Gain(lambda) = [EI_real(z=0) - EI_real(z=1)] - [EI_null(z=0) - EI_null(z=1)]

If Gain > 0 and exceeds bootstrap error, the topological evolution is
distinct from clustering evolution.

Graph construction: Fixed-Radius (percolation) graph with gravity weighting.
Units: All positions in cMpc/h (comoving). TNG300-1 box = 205/0.6774 cMpc/h.
"""

import os
import numpy as np
from scipy.spatial import cKDTree
from scipy.stats import entropy
from scipy.sparse import lil_matrix
import matplotlib.pyplot as plt

# --- CONFIGURATION ---
BOX_SIZE     = 205.0 / 0.6774       # TNG300-1 box in cMpc/h (~302.6)
N_SCALES     = 15
R_VALS       = np.linspace(2.0, 60.0, N_SCALES)
N_BOOTSTRAP  = 10
EPS          = 1e-12
SIGMA_THRESH = 2.0
DATA_DIR     = '/home/ubuntu/Dimensionality-Reduction-as-Cognitive-Fallacy'

# --- NORMALIZED EI OPERATOR ---
def get_normalized_ei(pos, mass, tree, r_link):
    """
    Compute normalized effective information at linking length r_link.
    Fixed-radius graph with asymmetric gravity weighting (m/d^2).
    Normalized by random-walk entropy H_rw.
    """
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
    """Compute EI_norm at each scale in R_VALS."""
    tree = cKDTree(pos, boxsize=BOX_SIZE)
    return np.array([get_normalized_ei(pos, mass, tree, r) for r in R_VALS])


def analyze_snapshot(label, pos, mass):
    """
    For one snapshot, compute:
      1. Real EI_norm curve
      2. Mass-shuffled null (mean + std over N_BOOTSTRAP)
      3. RGG null (Poisson with same density)
    """
    N = len(pos)
    print(f"  Analyzing {label} (N={N})...")

    # 1. Real topology
    print(f"    Computing real EI_norm ({N_SCALES} scales)...")
    ei_real = compute_ei_curve(pos, mass)

    # 2. Mass-shuffled null (geometry control)
    print(f"    Computing {N_BOOTSTRAP} mass-shuffled nulls...")
    ei_shuffles = []
    for b in range(N_BOOTSTRAP):
        m_rand = np.random.permutation(mass)
        ei_shuffles.append(compute_ei_curve(pos, m_rand))
        if (b + 1) % 5 == 0:
            print(f"      Bootstrap {b+1}/{N_BOOTSTRAP} done")
    ei_null_mean = np.mean(ei_shuffles, axis=0)
    ei_null_std  = np.std(ei_shuffles, axis=0)

    # 3. RGG null (clustering control)
    # Poisson-distributed points with same N, same box, same masses
    print(f"    Computing RGG (Poisson) null...")
    pos_rgg = np.random.rand(N, 3) * BOX_SIZE
    ei_rgg = compute_ei_curve(pos_rgg, mass)

    return {
        'ei_real':      ei_real,
        'ei_null_mean': ei_null_mean,
        'ei_null_std':  ei_null_std,
        'ei_rgg':       ei_rgg,
        'N':            N,
    }


def plot_double_difference(res_z0, res_z1, out_path):
    """
    Produce the two-panel figure:
      Left:  Absolute EI_norm evolution (real + RGG baseline)
      Right: Double Difference (Topological Gain)
    """
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))

    # --- Left: Absolute EI_norm ---
    ax = axes[0]
    ax.plot(R_VALS, res_z0['ei_real'], 'b-o', markersize=4, linewidth=2,
            label=f"z=0 Real (N={res_z0['N']})")
    ax.plot(R_VALS, res_z1['ei_real'], color='darkorange', marker='o',
            markersize=4, linewidth=2, label=f"z=1 Real (N={res_z1['N']})")
    ax.plot(R_VALS, res_z0['ei_rgg'], 'b:', alpha=0.5, linewidth=1.5,
            label='z=0 RGG Baseline')
    ax.plot(R_VALS, res_z1['ei_rgg'], color='darkorange', linestyle=':',
            alpha=0.5, linewidth=1.5, label='z=1 RGG Baseline')

    # Mass-shuffled null bands
    ax.fill_between(R_VALS,
                    res_z0['ei_null_mean'] - SIGMA_THRESH * res_z0['ei_null_std'],
                    res_z0['ei_null_mean'] + SIGMA_THRESH * res_z0['ei_null_std'],
                    color='blue', alpha=0.08, label=f'z=0 Null (+/-{SIGMA_THRESH:.0f}s)')
    ax.fill_between(R_VALS,
                    res_z1['ei_null_mean'] - SIGMA_THRESH * res_z1['ei_null_std'],
                    res_z1['ei_null_mean'] + SIGMA_THRESH * res_z1['ei_null_std'],
                    color='orange', alpha=0.08, label=f'z=1 Null (+/-{SIGMA_THRESH:.0f}s)')

    ax.set_xlabel('Linking Length (cMpc/h)', fontsize=12)
    ax.set_ylabel('EI_norm', fontsize=12)
    ax.set_title('Absolute EI_norm Evolution\nTNG300-1 (Top 1000 by Mass)', fontsize=12)
    ax.legend(fontsize=8, loc='best')
    ax.grid(True, alpha=0.3)

    # --- Right: Double Difference ---
    ax2 = axes[1]

    # Real gain
    delta_real = res_z0['ei_real'] - res_z1['ei_real']

    # Mass-shuffled null gain (mean)
    delta_null = res_z0['ei_null_mean'] - res_z1['ei_null_mean']

    # Bootstrap error on the null gain (propagated)
    delta_null_std = np.sqrt(res_z0['ei_null_std']**2 + res_z1['ei_null_std']**2)

    # Topological gain
    gain = delta_real - delta_null

    # RGG gain for reference
    delta_rgg = res_z0['ei_rgg'] - res_z1['ei_rgg']
    gain_vs_rgg = delta_real - delta_rgg

    ax2.plot(R_VALS, gain, 'g-o', linewidth=2.5, markersize=4,
             label='Topological Gain (vs Mass-Shuffled)')
    ax2.fill_between(R_VALS,
                     -SIGMA_THRESH * delta_null_std,
                     +SIGMA_THRESH * delta_null_std,
                     color='gray', alpha=0.2,
                     label=f'+/-{SIGMA_THRESH:.0f}s Null Uncertainty')
    ax2.plot(R_VALS, gain_vs_rgg, 'm--', linewidth=1.5, alpha=0.7,
             label='Topological Gain (vs RGG)')
    ax2.axhline(0, color='black', linestyle='--', alpha=0.8)

    # Check significance
    sig_mask = np.abs(gain) > SIGMA_THRESH * delta_null_std
    n_sig = int(np.sum(sig_mask))
    peak_gain = gain[np.argmax(np.abs(gain))]
    peak_scale = R_VALS[np.argmax(np.abs(gain))]

    ax2.set_xlabel('Linking Length (cMpc/h)', fontsize=12)
    ax2.set_ylabel('Gain (z=0 - z=1)', fontsize=12)
    ax2.set_title(f'Double Difference: Topological Gain\n'
                  f'Peak={peak_gain:.4f} at {peak_scale:.1f} cMpc/h, '
                  f'{n_sig}/{N_SCALES} scales >{SIGMA_THRESH:.0f}s',
                  fontsize=11)
    ax2.legend(fontsize=9, loc='best')
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    print(f"\nFigure saved to {out_path}")

    # Print summary
    print("\n=== DOUBLE DIFFERENCE RESULTS ===")
    print(f"  z=0 N: {res_z0['N']},  z=1 N: {res_z1['N']}")
    print(f"  Peak Topological Gain: {peak_gain:.6f} at {peak_scale:.1f} cMpc/h")
    print(f"  Scales exceeding {SIGMA_THRESH:.0f}s: {n_sig}/{N_SCALES}")
    print(f"  Mean |Gain|: {np.mean(np.abs(gain)):.6f}")
    print(f"  Mean |Gain vs RGG|: {np.mean(np.abs(gain_vs_rgg)):.6f}")

    # Per-scale detail
    print(f"\n  {'Scale':>8s}  {'Gain':>10s}  {'NullStd':>10s}  {'Sig':>5s}")
    for i, r in enumerate(R_VALS):
        sig_flag = "*" if sig_mask[i] else ""
        print(f"  {r:8.1f}  {gain[i]:10.6f}  {delta_null_std[i]:10.6f}  {sig_flag:>5s}")
    print("==================================")


if __name__ == "__main__":
    np.random.seed(42)

    # Load data
    pos_z0  = np.load(os.path.join(DATA_DIR, "tng300_z=0_pos.npy"))
    mass_z0 = np.load(os.path.join(DATA_DIR, "tng300_z=0_mass.npy"))
    pos_z1  = np.load(os.path.join(DATA_DIR, "tng300_z=1_pos.npy"))
    mass_z1 = np.load(os.path.join(DATA_DIR, "tng300_z=1_mass.npy"))

    print(f"Loaded z=0: {pos_z0.shape}, z=1: {pos_z1.shape}")
    print(f"Box size: {BOX_SIZE:.2f} cMpc/h")

    # Analyze both snapshots
    print("\n=== z=0 (Mature Web) ===")
    res_z0 = analyze_snapshot("z=0", pos_z0, mass_z0)

    print("\n=== z=1 (Early Web) ===")
    res_z1 = analyze_snapshot("z=1", pos_z1, mass_z1)

    # Plot
    out_path = os.path.join(DATA_DIR, 'double_difference.png')
    plot_double_difference(res_z0, res_z1, out_path)
