---

# Secular Evolution of Scale-Dependent Causal Structure in the Cosmic Web

**Zoverions**

*Independent Researcher in Complexity Theory*

---

## Abstract

Quantifying the emergence of macroscopic structure in cosmological simulations requires metrics that distinguish topological routing efficiency from simple spatial clustering. We introduce Scale-Dependent Effective Information (SDEI), a renormalized information-theoretic operator designed to measure the causal constraint of the cosmic web. We apply the SDEI operator to the 1,000 most massive subhalos in the IllustrisTNG300-1 simulation at $z=1$ and $z=0$. First, we validate the operator by showing that the cosmic web exhibits significantly lower effective information than density-matched random geometric graphs, confirming that gravitational collapse acts as a global constraint on information diffusion. Second, in a "Double Difference" analysis controlling for clustering evolution, we find no statistically significant topological gain at $N=1000$. The evolution of the causal signal is fully explained by the growth of the two-point correlation function. We establish this as a resolution floor: at the scale of massive clusters ($\sim 30$ cMpc mean separation), topological emergence is indistinguishable from density evolution.

## I. Introduction

Standard cosmological metrics (power spectrum, correlation function) characterize the clustering of matter. Network science offers tools to measure the routing architecture of these structures. Does the cosmic web evolve into a more efficient information-processing structure over time (topological modification), or does it merely get denser (clustering modification)?

To answer this, we employ a "Double Difference" protocol using the Scale-Dependent Effective Information (SDEI) operator. We normalize the predictive information of a diffusion process on the halo graph by the entropy of a random walk on the same degree sequence ($EI_{norm}$), isolating network topology from trivial degree-density effects.

## II. Methodology

We analyze the IllustrisTNG300-1 simulation, selecting the top 1,000 subhalos by mass at $z=1$ and $z=0$. This creates a sparse node distribution with a mean inter-halo separation of $\sim 30$ cMpc/h.

### SDEI Operator Definition

At a given linking length $\lambda$, we construct a **fixed-radius graph** where nodes are connected if their separation is less than $\lambda$. The graph is weighted by an asymmetric gravity term, where the edge from node $i$ to $j$ is $w_{ij} = m_j / d_{ij}^2$. The resulting adjacency matrix is made symmetric and then row-normalized to form a transition matrix $W$ for a diffusion process.

The operator is defined as:

$$
EI_{norm}(\lambda) = \frac{H_{eff} - H_{noise}}{H_{rw}}
$$

where $H_{eff}$ is the entropy of the stationary distribution of the diffusion process, $H_{noise}$ is the mean entropy of the transition probabilities from each node, and $H_{rw}$ is the entropy of a random walk on the same degree sequence, computed analytically from the node degrees.

### The Double Difference Protocol

To distinguish topological evolution from simple clustering growth, we define the Topological Gain:

$$
\text{Gain}(\lambda) = \Delta EI_{real} - \Delta EI_{null}
$$

where $\Delta EI = EI(z=0) - EI(z=1)$.

The Null: A mass-shuffled control that preserves the exact spatial positions of halos at each snapshot but randomizes their gravitational weights. This control captures the $EI$ growth driven purely by the steepening of the two-point correlation function.

## III. Results

### 1. Operator Validation: The Global Constraint of Gravity

Comparing the Real Web to a mass-preserving Random Geometric Graph (RGG) baseline reveals a massive deficit in effective information. The RGG produces high $EI_{norm}$ peaks ($>2.0$) at large scales due to uniform connectivity. The Real Web remains suppressed ($<0.3$). This confirms that gravity breaks the "small-world" randomness of a uniform field and validates that the SDEI operator is sensitive to spatial structure.

### 2. The Clustering Dominance ($N=1000$)

We observe a raw increase in $EI_{norm}$ at small scales ($2$ cMpc/h) from $z=1$ ($\sim 0.15$) to $z=0$ ($\sim 0.27$). However, the Double Difference analysis returns a null result. The mass-shuffled null reproduces this growth trend with high fidelity. The Topological Gain is $\approx -0.04 \pm \text{noise}$, effectively zero. The observed increase in causal constraint is fully explained by the evolution of the underlying spatial distribution (clustering). At the resolution of massive clusters, the specific assignment of masses to nodes does not create a distinct topological signal beyond what is predicted by the density field itself.

## IV. Discussion: The Resolution Floor

This null result establishes a critical boundary condition for topological cosmology. At $N=1000$ in a 300 Mpc volume, the graph saturates (becomes fully connected) at linking lengths $> 6$ cMpc/h. The "Sponge Regime"—where filamentary routing occurs—typically exists between $10-30$ Mpc. However, because the mean separation of clusters ($\sim 30$ Mpc) is comparable to this scale, the filamentary structure is undersampled and washed out by all-to-all connectivity.

Our spatial shuffle diagnostic confirms the operator is sensitive to spatial structure, but the Double Difference null result indicates that at this resolution, the mass-position correlations are not yet strong enough to be distinguished from the underlying density field. The saturation of the graph at small linking lengths is consistent with a resolution artifact, but we cannot rule out that the null is genuine. To detect the topological signal of the cosmic web, future analysis must resolve the filamentary bridges themselves, not just the nodes. We predict that increasing the sample size to $N \geq 15,000$ (reducing mean separation to $\sim 9$ cMpc/h) is the necessary threshold to disentangle topology from density. Until then, the hypothesis of topological evolution remains distinct from, but subordinate to, density evolution.

## V. Data Availability

All scripts, catalogs, and analysis pipelines are available in the repository at [github.com/Zoverions/Secular-Evolution-of-Causal-Structure](https://github.com/Zoverions/Secular-Evolution-of-Causal-Structure).

## References

[1] Nelson, D., et al. (2019). The IllustrisTNG simulations. *Computational Astrophysics and Cosmology*.

[2] Hoel, E. P. (2013). Quantifying causal emergence. *PNAS*.

[3] Springel, V., et al. (2018). First results from the IllustrisTNG simulations. *MNRAS*.
