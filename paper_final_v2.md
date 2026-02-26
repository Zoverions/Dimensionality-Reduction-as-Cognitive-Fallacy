# Secular Evolution of Scale-Dependent Causal Structure: A Resolution Sensitivity Analysis of IllustrisTNG

**Zoverions**

*Independent Researcher in Complexity Theory*

### Abstract

Quantifying the emergence of macroscopic structure in cosmological simulations requires metrics that distinguish topological routing efficiency from simple spatial clustering. We introduce **Scale-Dependent Effective Information (SDEI)**, a renormalized information-theoretic operator designed to measure the causal constraint of the cosmic web. We report two primary findings. First, the cosmic web exhibits significantly lower effective information than density-matched random geometric graphs, confirming that gravitational collapse imposes genuine causal constraint. Second, in a resolution-ladder analysis, the double-difference topological gain is statistically consistent with zero at N=1,000 (mean separation ∼30 cMpc/h), but becomes significantly negative (−0.12 ± 0.02) at N=5,000 (∼18 cMpc/h). This negative signal—robust across λ ≳ 10 cMpc/h—reveals that filamentary topology suppresses the growth of causal routing efficiency beyond what clustering evolution alone predicts. We interpret this as emerging topological rigidity in the cosmic web. The result is resolution-dependent: early-web (z=1) EI_norm at large scales rises markedly with N, indicating that coarse sampling washes out filamentary structure. We conclude that topological emergence in the cosmic web is masked at cluster-only resolution and becomes detectable once the filamentary network is adequately sampled (N ≳ 5,000). Scripts and catalogs are provided for extension to N=15k on HPC.

### I. Introduction

Standard cosmological metrics (power spectrum, correlation function) characterize the *clustering* of matter. Network science offers tools to measure the *routing architecture* of these structures. Does the cosmic web evolve into a more efficient information-processing structure over time (topological modification), or does it merely get denser (clustering modification)?

To answer this, we employ a "Double Difference" protocol using the **Scale-Dependent Effective Information (SDEI)** operator. We normalize the predictive information of a diffusion process on the halo graph by the entropy of a random walk on the same degree sequence ($EI_{norm}$), isolating network topology from trivial degree-density effects.

### II. Methodology

#### 1. Graph Construction and Operator Definition

We analyze the **IllustrisTNG300-1** simulation, selecting the top $N$ subhalos by mass at $z=1$ and $z=0$.

For a set of $N$ halos with positions $\mathbf{x}$ and masses $M$, we construct a **fixed-radius spatial graph** $G_\lambda$ at linking length $\lambda$. The adjacency matrix $A$ is defined with **asymmetric gravity-weighted edges**:

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

![Resolution Ladder](resolution_ladder.png)

**Figure 1: Resolution Ladder for Scale-Dependent Effective Information in TNG300-1 (N=1k vs N=5k).** *Left:* Double-difference topological gain (real minus mass-shuffled null; signal above ±2σ null band). Blue: N=1,000 (mean separation 30 cMpc/h); orange: N=5,000 (18 cMpc/h). *Center:* Absolute EI_norm at z=0 (snapshot 99, mature web). *Right:* Absolute EI_norm at z=1 (snapshot 50, early web). The annotation highlights the resolution effect at large linking lengths. Note: N=15k (mean separation ∼13 cMpc/h) requires an HPC cluster; scripts provided in the repository.

#### 1. Validation: The Global Constraint of Gravity

Comparing the Real Web to a mass-preserving Random Geometric Graph (RGG) baseline reveals a massive deficit in effective information.

* **Observation:** The RGG produces high $EI_{norm}$ peaks ($>1.5$) at large scales due to uniform connectivity. The Real Web remains suppressed ($<0.3$).
* **Interpretation:** This result validates the operator’s sensitivity. Gravity suppresses the near-random routing efficiency of a uniform field. By clustering matter into isolated potentials, it constrains the diffusion of information relative to the random baseline.

#### 2. The Clustering Dominance ($N=1000$)

We observe a raw increase in $EI_{norm}$ at small scales ($2$ cMpc/h) from $z=1$ ($\sim 0.15$) to $z=0$ ($\sim 0.27$). However, the **Double Difference analysis returns a null result.**

* **Finding:** The mass-shuffled null reproduces this growth trend with high fidelity. The Topological Gain is $\approx -0.04 \pm \text{noise}$, effectively zero.
* **Conclusion:** The observed increase in causal constraint is fully explained by the evolution of the underlying spatial distribution (clustering). At the resolution of massive clusters, the specific assignment of masses to nodes does not create a distinct topological signal beyond what is predicted by the density field itself.

#### 3. Resolution Sensitivity: Negative Topological Gain Emerges at N=5,000

When the halo sample is increased to the 5,000 most massive subhalos, the double-difference signal departs markedly from zero (Fig. 1, left, orange). The topological gain stabilizes at ≈ −0.12 across 10–60 cMpc/h, statistically significant above the null band. Concurrently, EI_norm(z=1) at large scales rises by ∼80 % relative to the N=1k case (Fig. 1, right), while the mature-web values (z=0) remain insensitive to resolution.

This demonstrates that the N=1k null result is a resolution artifact: coarse sampling (mean inter-halo distance ∼30 cMpc/h) saturates the graph before the filamentary “sponge regime” (∼10–30 cMpc/h) is properly resolved. At N=5k the filamentary bridges become visible, and the real gravitational weighting produces a *net suppression* of causal-information growth relative to clustering alone—evidence of topological rigidity.

### IV. Discussion

At N=1,000 the analysis returns a clean null, showing that topological emergence is indistinguishable from density evolution at cluster scales. The resolution ladder, however, reveals that this null is resolution-limited. At N=5,000 the double-difference becomes robustly negative, while early-web EI_norm at large linking lengths increases with sample size. Together these findings indicate that the cosmic web self-organizes into filamentary structures that impose additional causal constraint beyond two-point clustering. The negative gain implies an optimization: mass assignment along filaments channels information more efficiently (lower normalized diffusive entropy) than a random shuffling of the same density field would allow.

We cannot yet claim that the cosmic web’s causal structure is *entirely* reducible to its two-point statistics; rather, the statistics alone miss the directed-routing architecture that only appears once filaments are resolved. Extension to N≥15k (mean separation ∼13 cMpc/h) will test whether the negative signal strengthens or saturates, and whether a positive topological gain eventually emerges at even finer scales.

### V. Data Availability

All scripts, catalogs, and analysis pipelines are available in the repository `Zoverions/Secular-Evolution-Causal-Structure`.

---

**References**

1. Hoel, E. P., Albantakis, L., & Tononi, G. (2013). Quantifying causal emergence shows that macro can beat micro. *PNAS*, 110(49), 19790-19795.
2. Nelson, D., et al. (2019). The IllustrisTNG simulations: public data release. *Computational Astrophysics and Cosmology*, 6(1), 1-19.
3. Pearl, J. (2009). *Causality: Models, Reasoning, and Inference*. Cambridge University Press.
4. Springel, V., et al. (2018). First results from the IllustrisTNG simulations: matter and galaxy clustering. *MNRAS*, 475(5), 676-698.
