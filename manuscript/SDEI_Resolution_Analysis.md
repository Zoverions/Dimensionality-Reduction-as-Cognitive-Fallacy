# Secular Evolution of Scale-Dependent Causal Structure: A Resolution Sensitivity Analysis of IllustrisTNG

**Zoverions**
*Independent Researcher in Complexity Theory*

### Abstract
Quantifying the emergence of macroscopic structure in cosmological simulations requires metrics that distinguish topological routing efficiency from simple spatial clustering. We introduce **Scale-Dependent Effective Information (SDEI)**, a renormalized information-theoretic operator designed to measure the causal constraint of the cosmic web. In this pilot study, we apply the SDEI operator to the 1,000 most massive subhalos in the **IllustrisTNG300-1** simulation at $z=1$ and $z=0$.

We report two primary findings. First, we validate the operator by demonstrating that the cosmic web exhibits significantly *lower* effective information than density-matched random geometric graphs, confirming that the metric correctly identifies gravitational collapse as a constraint on information diffusion relative to a random field. Second, in a "Double Difference" analysis controlling for clustering evolution, we find no statistically significant topological gain at $N=1000$. The evolution of the causal signal is fully explained by the growth of the two-point correlation function. We conclude that at the coarse resolution of massive clusters (~30 cMpc mean separation), topological emergence is indistinguishable from density evolution. We discuss whether this represents a genuine null result or a resolution artifact due to graph saturation, and propose a specific threshold ($N \ge 5,000$) where topological rigidity is predicted to emerge.

### I. Introduction
Standard cosmological metrics (power spectrum, correlation function) characterize the *clustering* of matter. Network science offers tools to measure the *routing architecture* of these structures. Does the cosmic web evolve into a more efficient information-processing structure over time (topological modification), or does it merely get denser (clustering modification)?

To answer this, we employ a "Double Difference" protocol using the **Scale-Dependent Effective Information (SDEI)** operator. We normalize the predictive information of a diffusion process on the halo graph by the entropy of a random walk on the same degree sequence ($EI_{norm}$), isolating network topology from trivial degree-density effects.

### II. Methodology

#### 1. Graph Construction and Operator Definition
We analyze the **IllustrisTNG300-1** simulation, selecting the top $N$ subhalos by mass at $z=1$ and $z=0$.

For a set of $N$ halos with positions $\mathbf{x}$ and masses $M$, we construct a **fixed-radius spatial graph** $G_\lambda$ at linking length $\lambda$. The adjacency matrix $A$ is defined with asymmetric gravity-weighted edges:

$$A_{ij} = \begin{cases} \frac{M_j}{|\mathbf{x}_i - \mathbf{x}_j|^2} & \text{if } |\mathbf{x}_i - \mathbf{x}_j| < \lambda \\ 0 & \text{otherwise} \end{cases}$$

This graph is converted into a row-stochastic Markov transition matrix $W$ where $W_{ij} = A_{ij} / \sum_k A_{ik}$.

**Effective Information ($EI$)** is calculated as the mutual information of the one-step transition process on $W$ under a uniform intervention:

$$EI(W) = H(\langle W_i \rangle) - \langle H(W_i) \rangle$$

where $H(\cdot)$ is the Shannon entropy, $\langle W_i \rangle$ is the average transition vector (the effect), and $\langle H(W_i) \rangle$ is the average row entropy (the noise).

**Normalization:** To control for the natural increase in diffusion entropy caused by higher edge density, we define the **Normalized Effective Information**:

$$EI_{norm}(\lambda) = \frac{EI(W)}{H_{rw}(W)}$$

where $H_{rw}$ is the entropy of the stationary distribution $\pi$ of a random walk on the graph ($\pi W = \pi$), computed via power iteration to account for the directed nature of the gravity-weighted graph.

#### 2. The Double-Difference Protocol
To rigorously distinguish topological evolution from simple clustering growth, we define the **Topological Gain**:

$$\text{Gain}(\lambda) = \Delta EI_{real} - \Delta EI_{null}$$

where $\Delta EI = EI(z=0) - EI(z=1)$.
* **The Null:** A mass-shuffled control that preserves the exact spatial positions of halos at each snapshot but randomizes their gravitational weights. This control captures the $EI$ growth driven purely by the steepening of the two-point correlation function. We note that the samples at $z=1$ and $z=0$ represent the most massive objects at each epoch rather than a fixed lineage; the mass-shuffled null handles this population difference internally.

### III. Results

#### 1. Validation: The Global Constraint of Gravity
Comparing the Real Web to a mass-preserving Random Geometric Graph (RGG) baseline reveals a massive deficit in effective information.
* **Observation:** The RGG produces high $EI_{norm}$ peaks (>1.5) at large scales due to uniform connectivity. The Real Web remains suppressed (<0.3).
* **Interpretation:** This result validates the operator's sensitivity. Gravity suppresses the near-random routing efficiency of a uniform field. By clustering matter into isolated potentials, it constrains the diffusion of information relative to the random baseline.

#### 2. The Clustering Dominance ($N=1000$)
We observe a raw increase in $EI_{norm}$ at small scales (2 cMpc/h) from $z=1$ (~0.15) to $z=0$ (~0.27). However, the **Double Difference analysis returns a null result.**
* **Finding:** The mass-shuffled null reproduces this growth trend with high fidelity. The Topological Gain is $\approx -0.04 \pm \text{noise}$, effectively zero.
* **Conclusion:** The observed increase in causal constraint is fully explained by the evolution of the underlying spatial distribution (clustering). At the resolution of massive clusters, the specific assignment of masses to nodes does not create a distinct topological signal beyond what is predicted by the density field itself.

### IV. Discussion
This study establishes a critical boundary condition for topological cosmology. At $N=1000$ in a 300 Mpc volume, the graph saturates (becomes fully connected) at linking lengths > 6 cMpc/h. The "Sponge Regime"—where filamentary routing occurs—typically exists between 10–30 Mpc. However, because the mean separation of clusters (~30 Mpc) is comparable to this scale, the filamentary structure is undersampled and washed out by all-to-all connectivity.

We cannot yet rule out that the null result is genuine—that is, that the cosmic web's causal structure is entirely reducible to its two-point clustering statistics. However, the saturation of the graph at small scales suggests that this is a resolution artifact. To definitively distinguish between a genuine null and a resolution limit, future work must extend this analysis to $N \ge 5,000$ (reducing mean separation to ~18 cMpc/h) to resolve the filamentary bridges connecting the nodes.

### V. Data Availability
All scripts, catalogs, and analysis pipelines are available in the repository `Zoverions/Secular-Evolution-Causal-Structure`.

***

**References**
1. Hoel, E. P., Albantakis, L., & Tononi, G. (2013). Quantifying causal emergence shows that macro can beat micro. *PNAS*, 110(49), 19790-19795.
2. Nelson, D., et al. (2019). The IllustrisTNG simulations: public data release. *Computational Astrophysics and Cosmology*, 6(1), 1-19.
3. Pearl, J. (2009). *Causality: Models, Reasoning, and Inference*. Cambridge University Press.
4. Springel, V., et al. (2018). First results from the IllustrisTNG simulations: matter and galaxy clustering. *MNRAS*, 475(5), 676-698.