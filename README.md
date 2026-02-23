# Causal Resurgence Toy Model

This repository contains the Python implementation of the "Causal Resurgence" research proposal: **"Dimensionality Reduction as Cognitive Fallacy: A Topological Hypothesis and Research Proposal"**.

The goal of this project is to simulate and measure "Causal Resurgence" — a phenomenon where macroscopic network topology exerts top-down causal power over its microscopic substrate. This is operationalized using the **Causal Renormalization Group (CRG)** framework.

## Introduction

The central hypothesis is that certain complex systems (like biological networks or the cosmic web) exhibit **Causal Resurgence ($\beta_C > 0$)**. This occurs when the effective causal power of a system *increases* as we coarse-grain to larger scales, rather than decaying as reductionism would predict.

This repository provides tools to test this hypothesis on synthetic topologies and, with an API key, on the IllustrisTNG cosmological simulation.

## Installation

The project requires Python 3 and the following scientific libraries:

```bash
pip install numpy scipy networkx matplotlib requests
```

## Usage

### 1. Synthetic Topology Ablation Study (Phase 1)
Run the ablation suite to compare three topological regimes: Random (Poisson), Clustered (Gaussian), and Filamentary (Cosmic Web Proxy).

```bash
python ablation_study.py
```

This generates `ablation_results.png`, a dual-panel visualization:
- **Left Panel**: Plots the Causal Beta Function ($\beta_C$) vs. Linking Length for all three topologies, each compared against its own mass-shuffled null hypothesis (±1 SD).
- **Right Panel**: A bar chart comparing the peak Causal Power of each topology, indicating statistical significance.

**Hypothesis**: The Filamentary topology should exhibit the strongest positive spike in $\beta_C$, indicating that the interconnected "web" structure is the primary driver of causal emergence.

### 2. TNG Temporal Evolution Pipeline (Phase 2)
If you have a valid IllustrisTNG API Key, you can run the temporal analysis to compare the "Early Web" ($z=1$) vs. the "Mature Web" ($z=0$).

1. Open `tng_pipeline.py`.
2. Insert your API Key in the `API_KEY` variable.
3. Run the script:

```bash
python tng_pipeline.py
```

This tests the "Smoking Gun" hypothesis: if Causal Resurgence is real, $\beta_C$ should increase as the universe evolves and structure becomes more defined.

## Methodology

### Effective Information (EI)
We measure the causal power of a macroscopic state $S_\lambda$ at scale $\lambda$ using Effective Information:

$$EI(\lambda) = I(Do(S_\lambda^t); S_\lambda^{t+1})$$

In this spatial graph context, we treat the graph topology at linking length $\lambda$ as the macroscopic state. We compute the mutual information of a random walker navigating the gravity-weighted network.

### Causal Beta Function ($\beta_C$)
The flow of causal power across scales is governed by $\beta_C$:

$$\beta_C(\lambda) \equiv \frac{d}{d \ln \lambda} EI(\lambda)$$

- **$\beta_C < 0$ (Reductionist Decay)**: Causal power is lost as we zoom out.
- **$\beta_C > 0$ (Causal Resurgence)**: Causal power *emerges* at the macroscopic scale.

## Results Interpretation

- **Filamentary vs. Random**: We expect the Filamentary topology to show a significant positive $\beta_C$ peak, while the Random topology should hover near zero or decay.
- **Null Hypothesis**: The shaded regions in the plots represent the behavior of the same spatial structure but with randomized masses. A signal is only considered "Emergent" if it significantly clears this noise floor.
