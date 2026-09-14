# From Distance to Nearness — Reproducibility Package

This repository contains the MATLAB code and data used to reproduce the numerical multifractal analysis in the manuscript:

**“From Distance to Nearness: Rethinking Geometry on the Tree of Life”**  
Kevin Hudnall

The numerical analysis implements the branching-process Random Iterated Function System (bRIFS) introduced by Hudnall & D'Souza (2025) and compares finite-depth simulations with the analytic multifractal spectrum derived in the manuscript.

## Repository structure

```text
.
├── README.md
├── LICENSE
├── code/
│   ├── BuildMultifractalTreeFn.m
│   ├── MakeRandomTree.m
│   ├── spectrum_from_Leaves.m
│   └── spectrum_from_S.m
└── data/
    ├── Scale_Matrix_25_Iter.mat
    └── Leaves_25_Iter.mat
