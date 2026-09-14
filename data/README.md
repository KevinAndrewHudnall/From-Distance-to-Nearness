# Data

This directory contains the saved numerical realization used for the multifractal analysis in:

**“From Distance to Nearness: Rethinking Geometry on the Tree of Life”**  
Kevin Hudnall

The saved data allow the manuscript figures to be reproduced without rerunning the full branching-process Random Iterated Function System (bRIFS) simulation.

## Files

### `Leaves_25_Iter.mat`

Contains the recursive lineage scales from the 25-iteration bRIFS realization used in the manuscript.

The MATLAB variable `Leaves` is a cell array in which

```matlab
Leaves{n}
```

contains the accumulated scale $B_v$ of every recursive lineage $v$ present at recursive depth $n$.

Thus:

- `Leaves{1}` contains the recursive lineages at depth 1;
- `Leaves{2}` contains the recursive lineages at depth 2;
- ...
- `Leaves{25}` contains the recursive lineages at depth 25.

Lineage multiplicity is preserved. If several distinct descendant lineages inherit the same scale, that scale appears multiple times in `Leaves{n}`.

This is required because the finite-depth multifractal calculation counts recursive lineages rather than distinct numerical scale values.

`Leaves_25_Iter.mat` is the primary data file used by `code/spectrum_from_Leaves.m`.

---

## Reproducing the manuscript figures

From the repository root, add the code directory to the MATLAB path and load the saved lineage data:

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

This produces two figures.

### Main figure

`Figure_1.pdf` contains:

1. a heatmap of the finite-depth contraction exponents across recursive depth; and
2. the depth-25 empirical multifractal spectrum compared with the analytic prediction.

For each recursive lineage $v$ at depth $n$,

$$
\alpha_n(v) = -\frac{\log B_v}{n}.
$$

If $N_n(\alpha,\Delta\alpha)$ is the number of depth-$n$ lineages in a bin centered at $\alpha$, the finite-depth spectrum is estimated by

$$
\widehat f_n(\alpha) = \frac{\log N_n(\alpha,\Delta\alpha)}{n\alpha}.
$$

For the uniform-contraction case used in the simulation,

$$
r\sim\mathrm{Uniform}(0,1),
$$

the analytic spectrum is

$$
f(\alpha)=\frac{\log(m\alpha)+1-\alpha}{\alpha},
$$

with

$$
m=\mathbb{E}[K]\approx1.657.
$$

### Convergence figure

`Figure_S1_Convergence.pdf` compares the finite-depth spectrum estimates at

$$
n=15,\quad 20,\quad 25
$$

with the analytic spectrum.

All three estimates are obtained from the same saved stochastic realization.

## Loading the scale matrix

The corresponding scale matrix can be loaded separately with

```matlab
load('data/Scale_Matrix_25_Iter.mat');
```

which loads the variable

```matlab
S
```

for the same bRIFS realization.

## Simulation parameters

The saved realization was generated using:

```text
Maximum offspring per internal branching event: 3
Maximum generations per finite branching tree: 2
Recursive bRIFS iterations:                    25
Contraction distribution:                      Uniform(0,1)
Expected recursive offspring count E[K]:       approximately 1.657
```

The contraction factors were generated using MATLAB's `rand`.

Because the model is stochastic, a newly generated realization will differ numerically from the data stored here. The saved files should therefore be used for exact reproduction of the manuscript figures.

## Important note on lineage multiplicity

Do not remove duplicate numerical values from `Leaves`.

Within a recursive replacement event, multiple terminal descendants can inherit the same newly realized scale. These descendants are nevertheless distinct recursive lineages.

For example, if a replacement event produces three terminal descendants at scale $B$, then all three occurrences of $B$ must remain in the next level of `Leaves`.

Deduplicating these values changes the lineage count
N_n(\alpha,\Delta\alpha)
$$

and therefore changes the estimated multifractal spectrum.
