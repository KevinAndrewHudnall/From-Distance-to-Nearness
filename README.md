# From Distance to Nearness — Reproducibility Package

This repository contains the MATLAB code and data used to reproduce the numerical multifractal analysis in the manuscript:

**“From Distance to Nearness: Rethinking Geometry on the Tree of Life”**  
Kevin Hudnall

The numerical example implements the branching-process Random Iterated Function System (bRIFS) introduced by Hudnall & D'Souza (2025) and compares finite-depth realizations of its recursive scale geometry with the analytic multifractal spectrum derived in the manuscript.

## Repository structure

```text
.
├── README.md
├── LICENSE
│
├── code/
│   ├── README.md
│   ├── BuildMultifractalTreeFn.m
│   ├── MakeRandomTree.m
│   └── spectrum_from_Leaves.m
│
└── data/
    ├── README.md
    └── Leaves_25_Iter.mat
```

## Overview

The bRIFS recursively couples finite branching realizations with inherited geometric contraction.

At each recursive replacement event:

1. a finite Galton–Watson branching tree is generated;
2. a random contraction factor is applied to the inherited scale;
3. each terminal leaf becomes a recursive lineage entering the next replacement step.

The numerical analysis follows the scale of every recursive lineage through this process.

If $B_v$ is the accumulated scale of a lineage $v$ at recursive depth $n$, its finite-depth contraction exponent is

$$
\alpha_n(v) =-\frac{\log B_v}{n}.
$$

Let

$$
N_n(\alpha,\Delta\alpha)
$$

denote the number of depth-$n$ lineages whose contraction exponents lie in a bin of width $\Delta\alpha$ centered at $\alpha$. The finite-depth multifractal spectrum is estimated directly by

$$
\widehat f_n(\alpha)=\frac{\log N_n(\alpha,\Delta\alpha)}{n\alpha}.
$$

For the minimal bRIFS used in the numerical example, the contraction factors satisfy

$$
r\sim\mathrm{Uniform}(0,1).
$$

Writing

$$
m=\mathbb{E}[K],
$$

where $K$ is the number of terminal descendants produced by one recursive replacement, the analytic spectrum is

$$
f(\alpha)=\frac{\log(m\alpha)+1-\alpha}{\alpha}.
$$

For the branching parameters used here,

$$
m\approx1.657.
$$

## Reproducing the manuscript figures

The exact recursive lineage scales used for the manuscript numerical example are provided in

```text
data/Leaves_25_Iter.mat
```

From the repository root, add the code directory to the MATLAB path and load the saved realization:

```matlab
addpath('code');

load('data/Leaves_25_Iter.mat');
```

Then run:

```matlab
m = 1.657;

[alphaCenters, fhatFinal, conv] = ...
    spectrum_from_Leaves( ...
        Leaves, ...
        m, ...
        'Figure_1.pdf', ...
        'Figure_S1_Convergence.pdf', ...
        30, ...
        [15 20 25]);
```

The function produces two figures.

### Main numerical figure

`Figure_1.pdf` contains:

1. a heatmap of the finite-depth contraction exponents across recursive depths; and
2. the depth-25 empirical multifractal spectrum compared with the analytic prediction.

### Finite-depth convergence figure

`Figure_S1_Convergence.pdf` compares the empirical spectrum at recursive depths

$$
n=10,\quad20,15,\quad20,\quad25
$$

with the analytic spectrum.

All displayed depths come from the same stochastic realization.

## Generating a new bRIFS realization

A new realization can be generated using:

```matlab
addpath('code');

MaxOffspring = 3;
MaxGens = 2;
ITERATIONS = 25;

[S, P, Leaves] = ...
    BuildMultifractalTreeFn( ...
        MaxOffspring, ...
        MaxGens, ...
        ITERATIONS);
```

The resulting `Leaves` object can then be analyzed with:

```matlab
m = 1.657;

[alphaCenters, fhatFinal, conv] = ...
    spectrum_from_Leaves( ...
        Leaves, ...
        m, ...
        'Figure_1.pdf', ...
        'Figure_S1_Convergence.pdf', ...
        30, ...
        [15 20 25]);
```

Because the bRIFS is stochastic, a newly generated realization will not reproduce the saved numerical realization exactly.

The supplied `data/Leaves_25_Iter.mat` file should therefore be used when exact reproduction of the manuscript figures is desired.

## Simulation parameters

The manuscript numerical illustration uses:

```text
Maximum offspring per internal branching event: 3
Maximum generations per finite branching tree: 2
Recursive bRIFS iterations:                    25
Contraction distribution:                      Uniform(0,1)
Expected recursive offspring count E[K]:       approximately 1.657
Spectrum bins:                                 30
Displayed convergence depths:                  15, 20, 25
```

The contraction factors are generated using MATLAB's `rand`.

The choice of 25 recursive iterations provides a compromise between finite-depth convergence and computational reproducibility. Because the number of recursive lineages grows exponentially with depth, substantially deeper realizations require rapidly increasing memory and computation.

## Code

The `code/` directory contains three MATLAB functions.

### `BuildMultifractalTreeFn.m`

Constructs a complete stochastic realization of the bRIFS.

It returns:

- `S` — reconstructed scale histories along complete lineage paths;
- `P` — reconstructed progeny histories along complete lineage paths;
- `Leaves` — level-wise recursive lineage scales with lineage multiplicity preserved.

The multifractal analysis in the manuscript uses `Leaves`.

### `MakeRandomTree.m`

Generates one finite Galton–Watson branching realization and returns the number of terminal leaves.

This terminal-leaf count is the recursive offspring variable $K$ used in the theoretical formulation.

### `spectrum_from_Leaves.m`

Computes the direct finite-depth multifractal spectrum from the recursive lineage scales and compares it with the analytic spectrum.

It produces both the main numerical figure and the finite-depth convergence figure.

See [`code/README.md`](code/README.md) for additional details.

## Data

The `data/` directory contains:

### `Leaves_25_Iter.mat`

The recursive lineage scales from the 25-iteration realization used for the manuscript numerical example.

`Leaves{n}` contains one accumulated scale value $B_v$ for every recursive lineage present at depth $n$.

Lineage multiplicity is preserved. If several distinct descendant lineages inherit the same scale, that value appears multiple times in `Leaves{n}`.

This multiplicity must be retained because the multifractal estimator counts lineages, not distinct numerical scale values.

See [`data/README.md`](data/README.md) for additional details.

## MATLAB requirements

The code was developed in MATLAB.

Generation of a new realization with `BuildMultifractalTreeFn.m` currently uses:

- `parfor`;
- `gpuArray`.

The full simulation may therefore require the MATLAB Parallel Computing Toolbox and a supported GPU.

The supplied `Leaves_25_Iter.mat` file allows the manuscript multifractal analysis to be reproduced without regenerating the full branching realization.

## Reference

Hudnall, K.  
**From Distance to Nearness: Rethinking Geometry on the Tree of Life.**  
Manuscript in preparation / under review.

Hudnall, K. and D'Souza, R. (2025).  
**What does the tree of life look like as it grows? Evolution and the multifractality of time.**  
*Journal of Theoretical Biology*, **607**, 112121.

## License

See [`LICENSE`](LICENSE).
