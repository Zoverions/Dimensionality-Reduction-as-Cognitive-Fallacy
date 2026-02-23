# @title CRG Ablation Study (Comparison of Topologies)
import numpy as np
from scipy.spatial import cKDTree
import matplotlib.pyplot as plt
from causal_resurgence_toy import get_effective_information

def generate_random_topology(N, L):
    """
    Pure Random (Poisson) Topology:
    Uniformly distributed points in the box.
    """
    return np.random.rand(N, 3) * L

def generate_clustered_topology(N, L, num_clusters=8, cluster_std=5.0):
    """
    Clustered (Gaussian Mixture) Topology:
    Points distributed around random cluster centers.
    """
    # Define cluster centers
    centers = np.random.rand(num_clusters, 3) * L

    # Assign points to clusters
    points = []
    points_per_cluster = N // num_clusters
    remainder = N % num_clusters

    for i in range(num_clusters):
        n_points = points_per_cluster + (1 if i < remainder else 0)
        # Gaussian distribution around center
        cluster_points = centers[i] + np.random.normal(0, cluster_std, (n_points, 3))
        # Handle periodic boundary conditions (wrap around)
        cluster_points = cluster_points % L
        points.append(cluster_points)

    return np.vstack(points)

def generate_filamentary_topology(N, L):
    """
    Filamentary Topology (Original Toy Model):
    Sinusoidal density perturbation.
    """
    pos = np.random.rand(N, 3) * L
    # Introduce structure: sinusoidal density perturbation
    pos[:, 0] = (pos[:, 0] + 5.0 * np.sin(4 * np.pi * pos[:, 1] / L)) % L
    return pos

def run_ablation_study():
    print("Running CRG Ablation Study...")

    # Parameters
    N, L = 800, 100.0
    r_vals = np.linspace(2.0, 25.0, 20)

    # Common Mass Distribution (Log-Normal)
    # Use the same mass distribution for all topologies to isolate geometric effects
    np.random.seed(42)
    mass_dist = np.random.lognormal(10, 1.5, N)

    topologies = {
        "Random (Poisson)": generate_random_topology(N, L),
        "Clustered (Gaussian Mixture)": generate_clustered_topology(N, L),
        "Filamentary (Cosmic Web Proxy)": generate_filamentary_topology(N, L)
    }

    results = {}

    plt.figure(figsize=(10, 6))

    colors = {
        "Random (Poisson)": "blue",
        "Clustered (Gaussian Mixture)": "orange",
        "Filamentary (Cosmic Web Proxy)": "green"
    }

    for name, pos in topologies.items():
        print(f"Processing: {name}...")

        # Build spatial tree
        tree = cKDTree(pos, boxsize=L)

        # Compute EI across scales
        ei_vals = []
        perc_vals = []
        for r in r_vals:
            ei, perc = get_effective_information(pos, mass_dist, tree, r, L)
            ei_vals.append(ei)
            perc_vals.append(perc)

        # Compute Beta_C
        beta_c = np.gradient(ei_vals, np.log(r_vals))
        results[name] = beta_c

        # Plot
        plt.plot(r_vals, beta_c, label=name, color=colors[name], linewidth=2)

        # Identify sponge window for Filamentary only (as reference) or for each?
        # Let's just plot the curves cleanly as requested.

    # Plot formatting
    plt.xlabel(r'Linking Length $\lambda$ [Mpc/h equivalent]')
    plt.ylabel(r'Causal Beta Function $\beta_C$')
    plt.title('Causal Resurgence: Topological Ablation Study')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.axhline(0, color='black', linewidth=0.5, linestyle='--') # Zero line for reference

    # Save figure
    plt.savefig('ablation_results.png', dpi=300)
    print("Ablation study completed. Figure saved to ablation_results.png")

if __name__ == "__main__":
    run_ablation_study()
