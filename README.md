# From Distance to Nearness — Reproducibility Package

This repository contains the data and minimal code needed to reproduce the numerical multifractal-spectrum calculation used in the manuscript:

**“From Distance to Nearness: Rethinking Geometry on the Tree of Life”**  
Kevin Hudnall

## Contents
- `BuildMultifractalTreeFn.*` — implements the branching-process Random Iterated Function System (bRIFS) of Hudnall & D'Souza 2025 to generate a multifractal phylogenetic tree.
- `MakeRandomTree.*` — builds a single finite Galton-Watson tree. Called by BuildMultifractalTreeFn to implement the bRIFS.
- `data/Scale_Matrix_25_Iter.*` — the scale matrix `S` used in the manuscript’s numerical illustration.
- `code/spectrum_from_S.*` — computes the multifractal spectrum `f(α)` from `S` (partition sums → τ(q) → Legendre transform).
- `spectrum_from_Leaves.*`

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
