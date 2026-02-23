# @title CRG Toy Model (No API Key Required)
import numpy as np
from scipy.spatial import cKDTree
from scipy.stats import entropy
from scipy.sparse import lil_matrix, csr_matrix
import networkx as nx
import matplotlib.pyplot as plt

def get_effective_information(pos, mass, tree, r_link, L):
    # Find pairs within r_link distance
    pairs = tree.query_pairs(r_link)
    N = len(pos)

    # Initialize adjacency matrix
    A = lil_matrix((N, N), dtype=np.float64)
    edges = []

    # Iterate over pairs to compute weights
    for i, j in pairs:
        # Calculate periodic distance vector
        delta = np.abs(pos[i] - pos[j])
        delta = np.where(delta > L/2, L - delta, delta)
        d = np.linalg.norm(delta)

        if d > 0:
            # Gravity-weighted edges
            A[i, j] = mass[j] / (d**2)
            A[j, i] = mass[i] / (d**2)
            edges.append((i, j))

    # Compute percolation fraction
    G = nx.Graph()
    G.add_nodes_from(range(N))
    G.add_edges_from(edges)

    if len(G.edges) > 0:
        largest_cc = max(nx.connected_components(G), key=len)
        perc_frac = len(largest_cc) / N
    else:
        perc_frac = 0.0

    if A.nnz == 0:
        return 0.0, perc_frac

    # Normalize matrix to create transition probability matrix W
    A = A.tocsr()
    degrees = np.array(A.sum(axis=1)).flatten()
    degrees[degrees == 0] = 1.0 # Avoid division by zero for isolated nodes
    W = A.multiply((1.0/degrees)[:, np.newaxis]).tocsr()

    # Compute row entropies (noise)
    # entropy(pk) computes Shannon entropy of distribution pk
    # We add 1e-12 to avoid log(0)
    row_ents = []
    for i in range(N):
        if W.indptr[i+1] > W.indptr[i]:
            row_data = W.data[W.indptr[i]:W.indptr[i+1]]
            row_ents.append(entropy(row_data + 1e-12))
        else:
            row_ents.append(0.0)

    H_noise = np.mean(row_ents)

    # Compute effective entropy of the mean output distribution
    mean_W = np.array(W.mean(axis=0)).flatten()
    H_eff = entropy(mean_W + 1e-12)

    # Effective Information = H_eff - H_noise
    return H_eff - H_noise, perc_frac

def run_toy_model():
    print("Running CRG Toy Model...")
    # Parameters
    N, L, N_BOOT = 800, 100.0, 20

    # Generate synthetic data
    np.random.seed(42)
    pos = np.random.rand(N, 3) * L
    # Introduce structure: sinusoidal density perturbation
    pos[:, 0] = (pos[:, 0] + 5.0 * np.sin(4 * np.pi * pos[:, 1] / L)) % L

    # Assign masses log-normally distributed
    mass_real = np.random.lognormal(10, 1.5, N)

    # Build spatial tree
    tree = cKDTree(pos, boxsize=L)

    # Define range of linking lengths (r)
    # Note: Appendix B says r_vals = np.linspace(2.0, 25.0, 20)
    r_vals = np.linspace(2.0, 25.0, 20)

    print(f"Computing EI for real topology across {len(r_vals)} scales...")
    results_real = [get_effective_information(pos, mass_real, tree, r, L) for r in r_vals]
    ei_real = [res[0] for res in results_real]
    perc_real = [res[1] for res in results_real]

    # Calculate Beta_C (derivative of EI with respect to log(r))
    # We use gradient
    beta_real = np.gradient(ei_real, np.log(r_vals))

    print(f"Running {N_BOOT} bootstrap controls...")
    beta_rands = []
    for b in range(N_BOOT):
        # Shuffle masses to destroy mass-topology correlation
        m_rand = np.copy(mass_real)
        np.random.shuffle(m_rand)

        # Compute EI for random control
        # Note: Position topology is same, only masses are shuffled
        results_rand = [get_effective_information(pos, m_rand, tree, r, L) for r in r_vals]
        ei_r = [res[0] for res in results_rand]

        # Calculate Beta_C for random control
        beta_rands.append(np.gradient(ei_r, np.log(r_vals)))

    # Compute mean and std dev of controls
    mu_null = np.mean(beta_rands, axis=0)
    std_null = np.std(beta_rands, axis=0)

    # --- Visualization Code ---
    print("Generating Figure 1...")
    plt.figure(figsize=(10, 6))

    # Identify "sponge topology" window roughly where percolation fraction is between 0.1 and 0.9
    # This is a heuristic based on typical percolation theory
    perc_indices = np.where((np.array(perc_real) > 0.1) & (np.array(perc_real) < 0.9))[0]
    if len(perc_indices) > 0:
        r_start = r_vals[perc_indices[0]]
        r_end = r_vals[perc_indices[-1]]
        plt.axvspan(r_start, r_end, color='gray', alpha=0.1, label='Sponge Topology Window')

    # Plot Null Hypothesis (Control)
    plt.plot(r_vals, mu_null, color='red', label='Null Hypothesis (Random Mass)', linestyle='--')
    plt.fill_between(r_vals, mu_null - std_null, mu_null + std_null, color='red', alpha=0.2)

    # Plot Real Topology
    plt.plot(r_vals, beta_real, color='green', label='Real Topology (Gravity Clustered)', linewidth=2)

    # Labels and Title
    plt.xlabel(r'Linking Length $\lambda$ [Mpc/h equivalent]')
    plt.ylabel(r'Causal Beta Function $\beta_C$')
    plt.title('Causal Resurgence: Null Hypothesis Protocol')
    plt.legend()
    plt.grid(True, alpha=0.3)

    # Save figure
    plt.savefig('figure1.png', dpi=300)
    print("Figure saved to figure1.png")

if __name__ == "__main__":
    run_toy_model()
