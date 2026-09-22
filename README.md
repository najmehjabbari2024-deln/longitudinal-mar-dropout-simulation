# Longitudinal MAR Dropout Simulation

Monte Carlo simulation study comparing **LME, GEE, and GLMM** under monotone Missing-at-Random (MAR) dropout in longitudinal clinical trials.

This repository contains the R code and selected results associated with the following study:

**Comparing Linear Mixed Models, Generalized Estimating Equations, and Generalized Linear Mixed Models Under Monotone Missing-at-Random Dropout in Longitudinal Clinical Trials: A Monte Carlo Simulation Study**

**Preprint:** 10.20944/preprints202609.1550.v1

## Structure

```text
notebooks/    R scripts for data generation, dropout, model fitting,
              simulation, metrics, and output generation

results/      Selected simulation results, figures, and summary tables
```

## Reproducibility

The main pipeline is:

```text
01 → 02 → 03 → 04 → 05 → 06
```

and can be run using:

```r
source("notebooks/run_pipeline.R")
```

For the full study design, methodology, and results, please see the associated paper.
