# Simulation studies

This folder contains code and results for the simulation studies in the Supplement of [*Local graph estimation: Interpretable network discovery for complex data*](https://github.com/omelikechi/localgraph-paper). The simulations evaluate the performance of PFS and benchmark it against existing graph estimation methods under both linear and nonlinear data-generating processes.

---

## Contents

- `simulation.py`: Runs the main simulation studies from the Supplement (100 trials each)
- `analyze_simulation_results.py`: Loads and summarizes simulation output into tables of TPR and FDR by radius
- `illustration.py`: Generates **Figures 1a-1d** in the main text (visual comparison of PFS vs baseline methods in a single trial)
- `simulate_block.py`: Generates synthetic data with block-structured precision matrices
- `methods.py`: Wrapper functions for running baseline methods (graphical lasso, nodewise lasso, etc.)
- `results/`: Contains results from simulation runs

---

## Simulation studies

### Linear setting

- Each trial draws $n=100$ samples from a $p=200$-dimensional multivariate normal distribution.
- The precision matrix has three blocks representing:
  1. the target variable $X_1$,
  2. its four true neighbors (diagonal subblock),
  3. and the remaining variables.
- True neighbors of $X_1$ are connected to the rest of the graph, but not to each other.
- Edge weights are randomly sampled and the matrix is scaled to have eigenvalues between 0.01 and 10.

### Nonlinear setting

- Same setup as the linear case, but $X_1$ is redefined as a nonlinear function of its neighbors:
  $$
  X_1 = \sum_{j=2}^5 \exp(-X_j^2 / 2) + \varepsilon
  $$
  where $\varepsilon$ is Gaussian noise scaled to achieve a signal-to-noise ratio (SNR) of 4.
- This mimics scenarios where both high and low covariate levels are harmful (e.g., in medicine or epidemiology).

---

## Included methods

- **PFS**: Pathwise Feature Selection (with IPSS as the base selector)
- **GLasso**: Graphical lasso (via `huge` R package)
- **NLasso**: Nodewise lasso (via `huge`)
- **BNWSL**, **DSGL**, **DSNWSL**, **GFCL**, **GFCSL**: Additional baseline methods from the `silggm` R package

---

## Results

- `results/simulation_results_linear.pkl`: Output for the linear setting (100 trials)
- `results/simulation_results_nonlinear.pkl`: Output for the nonlinear setting (100 trials)
- Results include edge recovery across radii 1–4, allowing comparison of true and false positives

---

## Figure illustration

- `illustration.py` produces **Figures 1a-1d** from the main text, comparing PFS, graphical lasso, and nodewise lasso on a single simulated example.
- Results are deterministic on a single machine using `random_seed = 168`, but may vary across platforms.

---

## Notes

- The `simulate_block.py` script provides the generative model used in both simulation studies and for Figure 1.
- Simulation results match the tables reported in the Supplement (see Tables S3–S6).
