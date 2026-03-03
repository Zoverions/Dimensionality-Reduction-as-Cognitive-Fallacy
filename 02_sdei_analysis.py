"""
02_sdei_analysis.py
Core SDEI Operator and Resolution Ladder Analysis for TNG300-1.
Implements exact stationary distribution (power iteration) and uniform intervention.
"""

import numpy as np
from scipy.spatial import cKDTree
from scipy.stats import entropy
from scipy.sparse import lil_matrix, csr_matrix
import matplotlib.pyplot as plt
import os

# --- CONFIGURATION ---
R_VALS = np.linspace(2.0, 60.0, 20)
BOX_SIZE = 205.0 / 0.704 # TNG300-1 Box Size in cMpc/h
N_BOOTSTRAP = 10 # Increase to 50 for final paper figures

def get_stationary_distribution(W, tol=1e-9, max_iter=200):
    """Computes exact stationary distribution pi W = pi via Power Iteration."""
    N = W.shape[0]
    v = np.ones(N) / N # Uniform start
    for _ in range(max_iter):
        v_next = v @ W
        if np.linalg.norm(v_next - v) < tol:
            return v_next
        v = v_next
    return v

def get_exact_sdei(pos, mass, tree, r_link, box_size):
    """Computes Scale-Dependent Effective Information."""
    N = len(pos)
    pairs = tree.query_pairs(r_link)

    if len(pairs) == 0: return 0.0

    # 1. Construct Gravity-Weighted Graph (Asymmetric)
    A = lil_matrix((N, N), dtype=np.float64)
    for i, j in pairs:
        delta = np.abs(pos[i] - pos[j])
        delta = np.where(delta > box_size/2, box_size - delta, delta)
        d = np.linalg.norm(delta)
        if d > 0:
            # Gravity: F_ij ~ M_j / r^2 (Influence of j on i)
            A[i, j] = mass[j] / (d**2)
            A[j, i] = mass[i] / (d**2)

    A = A.tocsr()

    # 2. Row-Stochastic Transition Matrix W
    row_sums = np.array(A.sum(axis=1)).flatten()
    row_sums[row_sums == 0] = 1.0
    W = A.multiply((1.0 / row_sums)[:, np.newaxis]).tocsr()

    # 3. Effective Information (Uniform Intervention)
    W_mean = np.array(W.mean(axis=0)).flatten()
    H_eff = entropy(W_mean + 1e-12)

    row_entropies = []
    for i in range(N):
        ptr_start, ptr_end = W.indptr[i], W.indptr[i+1]
        if ptr_end > ptr_start:
            data = W.data[ptr_start:ptr_end]
            row_entropies.append(entropy(data + 1e-12))
        else:
            row_entropies.append(0.0)
    H_noise = np.mean(row_entropies)

    EI = H_eff - H_noise

    # 4. Normalization (H_rw) via exact stationarity
    pi = get_stationary_distribution(W)
    H_rw = entropy(pi + 1e-12)

    return EI / H_rw if H_rw > 0 else 0.0

def run_resolution_ladder(pos_z0, mass_z0, pos_z1, mass_z1):
    subsets = [1000, 5000, 15000]
    results = {}

    for n in subsets:
        print(f"\n--- Running Resolution Step N={n} ---")
        p0, m0 = pos_z0[:n], mass_z0[:n]
        p1, m1 = pos_z1[:n], mass_z1[:n]

        tree0 = cKDTree(p0, boxsize=BOX_SIZE)
        tree1 = cKDTree(p1, boxsize=BOX_SIZE)

        print("Calculating Real Web SDEI...")
        ei0 = [get_exact_sdei(p0, m0, tree0, r, BOX_SIZE) for r in R_VALS]
        ei1 = [get_exact_sdei(p1, m1, tree1, r, BOX_SIZE) for r in R_VALS]
        delta_real = np.array(ei0) - np.array(ei1)

        print("Calculating Mass-Shuffled Null...")
        null_deltas = []
        for b in range(N_BOOTSTRAP):
            m0_shuff = np.random.permutation(m0)
            m1_shuff = np.random.permutation(m1)
            e0_n = [get_exact_sdei(p0, m0_shuff, tree0, r, BOX_SIZE) for r in R_VALS]
            e1_n = [get_exact_sdei(p1, m1_shuff, tree1, r, BOX_SIZE) for r in R_VALS]
            null_deltas.append(np.array(e0_n) - np.array(e1_n))

        delta_null = np.mean(null_deltas, axis=0)
        gain = delta_real - delta_null

        results[n] = gain
        print(f"N={n} Peak Gain: {np.max(gain):.4f}")

    return results

if __name__ == "__main__":
    try:
        print("Loading TNG300 arrays...")
        p0 = np.load("data/tng300_z=0_pos.npy")
        m0 = np.load("data/tng300_z=0_mass.npy")
        p1 = np.load("data/tng300_z=1_pos.npy")
        m1 = np.load("data/tng300_z=1_mass.npy")

        ladder_results = run_resolution_ladder(p0, m0, p1, m1)

        plt.figure(figsize=(10, 6))
        for n, gain in ladder_results.items():
            plt.plot(R_VALS, gain, marker='o', label=f'N={n}')

        plt.axhline(0, color='k', linestyle='--', label='Null Hypothesis')
        plt.xlabel('Linking Length $\lambda$ (cMpc/h)')
        plt.ylabel('Topological Gain ($\Delta EI_{real} - \Delta EI_{null}$)')
        plt.title('SDEI Resolution Ladder: TNG300-1')
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.savefig("data/resolution_ladder_plot.png", dpi=300)
        print("Analysis complete. Plot saved to data/resolution_ladder_plot.png")

    except FileNotFoundError:
        print("Data files not found in data/. Please run 01_fetch_tng300.py first.")