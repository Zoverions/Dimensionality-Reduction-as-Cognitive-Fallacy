# @title CRG Toy Model Ablation Suite (Normalized, 2σ Threshold)
# Updated: falsification threshold upgraded to ±2σ per revised paper draft.
import numpy as np
from scipy.spatial import cKDTree
from scipy.stats import entropy
from scipy.sparse import lil_matrix, csr_matrix
import matplotlib.pyplot as plt

# --- CONFIGURATION ---
N_NODES = 1500
BOX_SIZE = 100.0
N_SCALES = 30
N_BOOTSTRAP = 15
R_MIN, R_MAX = 2.0, 25.0
EPS = 1e-12
SIGMA_THRESH = 2.0   # Updated: 2σ falsification threshold

# --- 1. TOPOLOGY GENERATORS ---
def generate_poisson(N, L):
    return np.random.rand(N, 3) * L

def generate_gaussian(N, L, n_clusters=5):
    centers = np.random.rand(n_clusters, 3) * L
    points = []
    per_cluster = N // n_clusters
    for center in centers:
        cluster_pts = np.random.normal(loc=center, scale=L/10, size=(per_cluster, 3))
        # Handle periodic boundary conditions roughly for cluster centers
        cluster_pts = cluster_pts % L
        points.append(cluster_pts)
    remainder = N - (per_cluster * n_clusters)
    if remainder > 0:
        points.append(np.random.rand(remainder, 3) * L)
    return np.vstack(points)

def generate_filamentary(N, L):
    pos = np.random.rand(N, 3) * L
    pos[:, 0] = (pos[:, 0] + (L/10) * np.sin(4 * np.pi * pos[:, 1] / L)) % L
    pos[:, 2] = (pos[:, 2] + (L/10) * np.cos(4 * np.pi * pos[:, 1] / L)) % L
    return pos

def assign_masses(N):
    return np.random.lognormal(mean=10, sigma=1.5, size=N)

# --- 2. NORMALIZED CORE OPERATORS ---
def compute_normalized_effective_information(pos, mass, r_link, L):
    N = len(pos)
    tree = cKDTree(pos, boxsize=L)
    pairs = tree.query_pairs(r_link, output_type='ndarray')
    if len(pairs) == 0: return 0.0

    A = lil_matrix((N, N), dtype=np.float64)
    i_idxs, j_idxs = pairs[:, 0], pairs[:, 1]

    delta = pos[i_idxs] - pos[j_idxs]
    delta = np.where(delta > L/2, L - delta, delta)
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

    # INFORMATION-THEORETIC NORMALIZATION
    node_probabilities = degrees / np.sum(degrees)
    H_rw = entropy(node_probabilities + EPS)

    return EI_raw / H_rw if H_rw > 0 else 0.0

def compute_normalized_beta_curve(pos, mass, L):
    r_vals = np.linspace(R_MIN, R_MAX, N_SCALES)
    ei_vals = [compute_normalized_effective_information(pos, mass, r, L) for r in r_vals]
    beta_vals = np.gradient(np.array(ei_vals), np.log(r_vals))
    return r_vals, np.array(ei_vals), beta_vals

# --- 3. ABLATION EXECUTION ---
def run_ablation_study():
    topologies = {
        "Poisson (Random)": generate_poisson,
        "Gaussian (Clusters)": lambda N, L: generate_gaussian(N, L),
        "Filamentary (Web)": generate_filamentary
    }
    results = {}
    print("Starting Normalized Ablation Study...")
    for name, gen_func in topologies.items():
        print(f"Processing {name}...")
        pos = gen_func(N_NODES, BOX_SIZE)
        mass = assign_masses(N_NODES)

        r, ei, beta = compute_normalized_beta_curve(pos, mass, BOX_SIZE)
        results[name] = {'r': r, 'beta': beta, 'ei': ei}

        beta_nulls = []
        for b in range(N_BOOTSTRAP):
            m_shuffle = np.random.permutation(mass)
            _, _, b_null = compute_normalized_beta_curve(pos, m_shuffle, BOX_SIZE)
            beta_nulls.append(b_null)

        results[f"{name} (Null)"] = {
            'r': r,
            'beta': np.mean(beta_nulls, axis=0),
            'std': np.std(beta_nulls, axis=0)
        }
    return results

# --- 4. VISUALIZATION ---
def plot_ablation(results):
    fig, ax = plt.subplots(1, 2, figsize=(16, 6))
    colors = {'Poisson': 'blue', 'Gaussian': 'orange', 'Filamentary': 'red'}
    base_names = ["Poisson (Random)", "Gaussian (Clusters)", "Filamentary (Web)"]

    for name in base_names:
        data = results[name]
        null_data = results[f"{name} (Null)"]
        color = colors[name.split()[0]]

        lower = null_data['beta'] - SIGMA_THRESH * null_data['std']
        upper = null_data['beta'] + SIGMA_THRESH * null_data['std']
        ax[0].fill_between(data['r'], lower, upper, color=color, alpha=0.15, label=f"{name} Null ±{SIGMA_THRESH:.0f}σ")
        ax[0].plot(data['r'], data['beta'], color=color, linewidth=2, label=f"{name} Real")

    ax[0].axhline(0, color='black', linestyle='--', alpha=0.5)
    ax[0].set_xlabel("Linking Length (r)")
    ax[0].set_ylabel("Normalized Causal Flow (Beta_C)")
    ax[0].set_title(f"Normalized Causal Resurgence Across Topologies (±{SIGMA_THRESH:.0f}σ)")
    ax[0].legend(loc='best', fontsize=8)
    ax[0].grid(True, alpha=0.3)

    peaks, names_short = [], []
    for name in base_names:
        peak_beta = np.max(results[name]['beta'])
        null_mean = results[f"{name} (Null)"]['beta']
        null_std = results[f"{name} (Null)"]['std']
        idx = np.argmax(results[name]['beta'])
        threshold = null_mean[idx] + SIGMA_THRESH * null_std[idx]
        significant = "Yes" if peak_beta > threshold else "No"
        peaks.append(peak_beta)
        names_short.append(f"{name.split()[0]}\n(Sig@{SIGMA_THRESH:.0f}σ: {significant})")

    ax[1].bar(names_short, peaks, color=['blue', 'orange', 'red'])
    ax[1].set_ylabel("Max Normalized Beta_C")
    ax[1].set_title(f"Peak Causal Power by Topology ({SIGMA_THRESH:.0f}σ Threshold)")
    ax[1].grid(True, alpha=0.3, axis='y')

    plt.tight_layout()
    plt.savefig('ablation_results.png')
    print("Ablation study completed. Figure saved to ablation_results.png")

if __name__ == "__main__":
    ablation_results = run_ablation_study()
    plot_ablation(ablation_results)
