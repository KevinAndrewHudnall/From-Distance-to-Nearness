# Multifractality in the Tree of Life — Reproducibility Package

This repository contains the data and minimal code needed to reproduce the numerical multifractal-spectrum calculation used in the manuscript:

**“Multifractality in the Tree of Life: A Branching-Process RIFS Proof”**  
Kevin Hudnall

## Contents
- `data/Scale_Matrix_20_Iter.*` — the scale matrix `S` used in the manuscript’s numerical illustration.
- `code/spectrum_from_S.*` — computes the multifractal spectrum `f(α)` from `S` (partition sums → τ(q) → Legendre transform).

## Generating the system (tree + scale matrix)
To generate new realizations of the branching-process RIFS and construct a scale matrix `S`, use the code in the companion repository:

- **the-living-tree-of-life** (primary codebase):  
  - function: `BuildMultifractalTreeFn`  
  - output: scale matrix `S`

(See that repository for installation and simulation settings.)

## Reproducing the spectrum from the provided S
1. Load `data/Scale_Matrix_20_Iter.*`
2. Run `code/spectrum_from_S.*`
3. The script outputs the estimated `τ(q)` and `f(α)` consistent with the manuscript.

## License
See `LICENSE`.
