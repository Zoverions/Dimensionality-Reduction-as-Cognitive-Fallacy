# Causal Resurgence Toy Model

This repository contains the Python implementation of the "Causal Resurgence Toy Model" as described in the research proposal **"Dimensionality Reduction as Cognitive Fallacy: A Topological Hypothesis and Research Proposal"**.

The goal of this project is to simulate and measure "Causal Resurgence" — a phenomenon where macroscopic network topology exerts top-down causal power over its microscopic substrate. This is operationalized using the **Causal Renormalization Group (CRG)** framework.

## Introduction

The central hypothesis is that certain complex systems (like biological networks or the cosmic web) exhibit **Causal Resurgence ($\beta_C > 0$)**. This occurs when the effective causal power of a system *increases* as we coarse-grain to larger scales, rather than decaying as reductionism would predict.

This toy model demonstrates this by comparing two scenarios:
1.  **Real Topology (Gravity-Clustered)**: A spatial graph where mass is distributed according to a structural pattern (sinusoidal density perturbation), creating a non-trivial causal topology.
2.  **Null Hypothesis (Random Mass)**: A control graph with the identical spatial positions but with masses randomly shuffled, destroying the correlation between mass and geometric structure.

If the "Real Topology" shows a positive spike in the Causal Beta Function ($\beta_C$) that significantly exceeds the Null Hypothesis, it indicates the emergence of a higher-order causal structure not present in the random control.

## Installation

The project requires Python 3 and the following scientific libraries:

```bash
pip install numpy scipy networkx matplotlib
```

## Usage

Run the simulation script directly:

```bash
python causal_resurgence_toy.py
```

This will:
1.  Generate a synthetic dataset of 800 nodes with a sinusoidal density structure.
2.  Compute the **Effective Information (EI)** across a range of spatial scales (linking lengths $\lambda$).
3.  Calculate the **Causal Beta Function ($\beta_C$)**, which measures the flow of causal power across scales.
4.  Run 20 bootstrap iterations of the Null Hypothesis (shuffled masses).
5.  Generate and save `figure1.png`, visualizing the results.

## Methodology

### 1. Effective Information (EI)
We measure the causal power of a macroscopic state $S_\lambda$ at scale $\lambda$ using Effective Information:

$$EI(\lambda) = I(Do(S_\lambda^t); S_\lambda^{t+1})$$

In this spatial graph context, we treat the graph topology at linking length $\lambda$ as the macroscopic state. We compute the mutual information of a random walker navigating the gravity-weighted network.

### 2. Causal Beta Function ($\beta_C$)
The flow of causal power across scales is governed by $\beta_C$:

$$\beta_C(\lambda) \equiv \frac{d}{d \ln \lambda} EI(\lambda)$$

- **$\beta_C < 0$ (Reductionist Decay)**: Causal power is lost as we zoom out (macro state is just noise).
- **$\beta_C \approx 0$ (Scale Invariance)**: Causal power is constant across scales.
- **$\beta_C > 0$ (Causal Resurgence)**: Causal power *emerges* at the macroscopic scale.

## Results

The script outputs `figure1.png`, which plots $\beta_C$ against the linking length $\lambda$.

- **Green Line**: The "Real Topology" (Gravity Clustered). A positive spike here indicates Causal Resurgence.
- **Red Line (with shaded band)**: The "Null Hypothesis" (Random Mass). This represents the baseline expectation for a system with no causal structure between mass and geometry.
- **Gray Band**: The "Sponge Topology Window", where the network is percolating but not fully connected, often where interesting topological phase transitions occur.

If the Green line significantly exceeds the Red band, we have successfully detected Causal Resurgence in the toy model.

## Peer Review & Ablation Study

Following a formal peer review, an ablation study was conducted to compare three distinct topological regimes. This experiment tests whether the "filamentary" structure of the cosmic web is uniquely capable of generating Causal Resurgence, compared to simpler geometries.

To run the ablation study:

```bash
python ablation_study.py
```

This generates `ablation_results.png`, comparing:

1.  **Random (Poisson) Topology**: Points uniformly distributed in space. Represents a structureless universe.
2.  **Clustered (Gaussian Mixture) Topology**: Points grouped in discrete blobs. Represents a universe with gravity but no large-scale filamentary network.
3.  **Filamentary (Cosmic Web Proxy)**: The sinusoidal structure used in the main toy model.

**Hypothesis**: The Filamentary topology should exhibit the strongest positive spike in $\beta_C$, indicating that the interconnected "web" structure is the primary driver of causal emergence, rather than simple clustering or random distribution.
