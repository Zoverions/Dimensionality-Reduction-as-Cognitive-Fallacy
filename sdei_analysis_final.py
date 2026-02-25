import numpy as np
from scipy.spatial import cKDTree
from scipy.stats import entropy
from scipy.sparse import lil_matrix, csr_matrix
import matplotlib.pyplot as plt
import time

# --- CONFIGURATION ---
R_VALS = np.linspace(2.0, 60.0, 20)
BOX_SIZE = 205.0 / 0.6774  # TNG300-1 (cMpc/h)
N_BOOTSTRAP = 10
EPS = 1e-12

def get_stationary_distribution(W, tol=1e-9, max_iter=100):
    """Computes exact stationary distribution pi W = pi via Power Iteration."""
    N = W.shape[0]
    v = np.ones(N) / N
    for _ in range(max_iter):
        v_next = v @ W
        s = v_next.sum()
        if s > 0:
            v_next /= s
        if np.linalg.norm(v_next - v) < tol:
            return v_next
        v = v_next
    return v

def get_exact_sdei(pos, mass, tree, r_link, box_size):
    """Computes Scale-Dependent Effective Information with exact stationarity."""
    N = len(pos)
    pairs = tree.query_pairs(r_link, output_type='ndarray')

    if len(pairs) == 0:
        return 0.0

    i_idxs, j_idxs = pairs[:, 0], pairs[:, 1]
    delta = pos[i_idxs] - pos[j_idxs]
    delta = np.where(delta > box_size / 2, delta - box_size, delta)
    delta = np.where(delta < -box_size / 2, delta + box_size, delta)
    dists = np.linalg.norm(delta, axis=1)

    mask = dists > 0
    i_idxs, j_idxs, dists = i_idxs[mask], j_idxs[mask], dists[mask]

    if len(i_idxs) == 0:
        return 0.0

    w_ij = mass[j_idxs] / (dists**2)
    w_ji = mass[i_idxs] / (dists**2)

    A = lil_matrix((N, N), dtype=np.float64)
    for k in range(len(i_idxs)):
        A[i_idxs[k], j_idxs[k]] = w_ij[k]
        A[j_idxs[k], i_idxs[k]] = w_ji[k]
    A = A.tocsr()

    row_sums = np.array(A.sum(axis=1)).flatten()
    row_sums[row_sums == 0] = 1.0
    W = A.multiply((1.0 / row_sums)[:, np.newaxis]).tocsr()

    W_mean = np.array(W.mean(axis=0)).flatten()
    H_eff = entropy(W_mean + EPS)

    row_entropies = np.zeros(N)
    for i in range(N):
        start, end = W.indptr[i], W.indptr[i+1]
        if end > start:
            row_entropies[i] = entropy(W.data[start:end] + EPS)
    H_noise = np.mean(row_entropies)
    EI = H_eff - H_noise

    pi = get_stationary_distribution(W)
    H_rw = entropy(pi + EPS)

    return EI / H_rw if H_rw > 0 else 0.0

def run_analysis(pos0, mass0, pos1, mass1):
    tree0 = cKDTree(pos0, boxsize=BOX_SIZE)
    tree1 = cKDTree(pos1, boxsize=BOX_SIZE)

    print("  Computing Real Evolution...")
    ei0 = [get_exact_sdei(pos0, mass0, tree0, r, BOX_SIZE) for r in R_VALS]
    ei1 = [get_exact_sdei(pos1, mass1, tree1, r, BOX_SIZE) for r in R_VALS]
    delta_real = np.array(ei0) - np.array(ei1)

    print("  Computing Null Evolution (Mass Shuffle)...")
    null_deltas = []
    for b in range(N_BOOTSTRAP):
        print(f"    Bootstrap {b+1}/{N_BOOTSTRAP}...")
        m0_s = np.random.permutation(mass0)
        m1_s = np.random.permutation(mass1)
        e0_n = [get_exact_sdei(pos0, m0_s, tree0, r, BOX_SIZE) for r in R_VALS]
        e1_n = [get_exact_sdei(pos1, m1_s, tree1, r, BOX_SIZE) for r in R_VALS]
        null_deltas.append(np.array(e0_n) - np.array(e1_n))

    delta_null_mean = np.mean(null_deltas, axis=0)
    delta_null_std = np.std(null_deltas, axis=0)
    gain = delta_real - delta_null_mean

    return gain, delta_null_std, ei0, ei1

if __name__ == "__main__":
    N_LADDER = [1000, 5000] #, 15000] # 15k requires HPC

    # Load full 15k datasets
    p0_full = np.load("tng300_z=0_pos_15000.npy")
    m0_full = np.load("tng300_z=0_mass_15000.npy")
    p1_full = np.load("tng300_z=1_pos_15000.npy")
    m1_full = np.load("tng300_z=1_mass_15000.npy")

    results = {}
    for n in N_LADDER:
        print(f"\n--- Running Analysis for N={n} ---")
        p0, m0 = p0_full[:n], m0_full[:n]
        p1, m1 = p1_full[:n], m1_full[:n]
        gain, std, ei0, ei1 = run_analysis(p0, m0, p1, m1)
        results[n] = {"gain": gain, "std": std, "ei0": ei0, "ei1": ei1}

    # Plotting
    fig, axes = plt.subplots(1, 3, figsize=(20, 6))
    colors = {1000: "#1f77b4", 5000: "#ff7f0e", 15000: "#2ca02c"}

    # Panel 1: Topological Gain
    for n in N_LADDER:
        res = results[n]
        axes[0].plot(R_VALS, res["gain"], "-o", color=colors[n], label=f"N={n}")
        axes[0].fill_between(R_VALS, -2 * res["std"], 2 * res["std"], color=colors[n], alpha=0.15)
    axes[0].axhline(0, color="k", linestyle="--")
    axes[0].set_title("Topological Gain (Double Difference)")
    axes[0].set_xlabel("Linking Length (cMpc/h)")
    axes[0].set_ylabel("Gain (Real - Null)")
    axes[0].legend()

    # Panel 2: EI_norm z=0
    for n in N_LADDER:
        axes[1].plot(R_VALS, results[n]["ei0"], "-o", color=colors[n], label=f"N={n}")
    axes[1].set_title("EI_norm at z=0")
    axes[1].set_xlabel("Linking Length (cMpc/h)")
    axes[1].set_ylabel("EI_norm")
    axes[1].legend()

    # Panel 3: EI_norm z=1
    for n in N_LADDER:
        axes[2].plot(R_VALS, results[n]["ei1"], "-o", color=colors[n], label=f"N={n}")
    axes[2].set_title("EI_norm at z=1")
    axes[2].set_xlabel("Linking Length (cMpc/h)")
    axes[2].set_ylabel("EI_norm")
    axes[2].legend()

    plt.tight_layout()
    plt.savefig("sdei_resolution_ladder.png", dpi=150)
    print("\nFigure saved as sdei_resolution_ladder.png")
