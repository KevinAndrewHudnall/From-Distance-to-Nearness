# Code

This directory contains the MATLAB functions used to generate the branching-process Random Iterated Function System (bRIFS) and reproduce the numerical multifractal analysis reported in:

**“From Distance to Nearness: Rethinking Geometry on the Tree of Life”**  
Kevin Hudnall

## Files

### `BuildMultifractalTreeFn.m`

Generates a realization of the bRIFS through recursive branching and scale contraction.

```matlab
[S, P, Leaves] = BuildMultifractalTreeFn(MaxOffspring, MaxGens, ITERATIONS);
```

Inputs:

- `MaxOffspring` — maximum offspring count for an internal branching event.
- `MaxGens` — maximum number of generations in each finite branching realization.
- `ITERATIONS` — number of recursive bRIFS replacement steps.

Outputs:

- `S` — scale matrix reconstructed along complete lineage paths.
- `P` — progeny-count matrix reconstructed along complete lineage paths.
- `Leaves` — cell array containing the scale of every recursive lineage present at each recursive depth.

`Leaves` preserves lineage multiplicity. Distinct descendant lineages therefore remain separate even when they inherit identical scale values. This is essential for the direct finite-depth multifractal estimator used in the manuscript.

For the numerical example in the manuscript:

```matlab
MaxOffspring = 3;
MaxGens = 2;
ITERATIONS = 25;

[S, P, Leaves] = ...
    BuildMultifractalTreeFn(MaxOffspring, MaxGens, ITERATIONS);
```

---

### `MakeRandomTree.m`

Generates one finite Galton–Watson branching realization.

```matlab
K = MakeRandomTree(MaxOffspring, MaxGens);
```

The returned value is the number of terminal leaves in the finite tree. In the notation of the manuscript, this is the recursive offspring count \(K\).

For the manuscript settings:

```matlab
MaxOffspring = 3;
MaxGens = 2;
```

the expected recursive offspring count is approximately

\[
m=\mathbb{E}[K]\approx 1.657.
\]

`BuildMultifractalTreeFn.m` calls `MakeRandomTree.m` once for each recursive lineage replacement.

---

### `spectrum_from_Leaves.m`

Performs the numerical multifractal analysis used in the manuscript.

For every recursive lineage \(v\) at depth \(n\), the function computes

\[
\alpha_n(v)
=
-\frac{\log B_v}{n}.
\]

If \(N_n(\alpha,\Delta\alpha)\) is the number of depth-\(n\) lineages in a bin centered at \(\alpha\), the finite-depth spectrum is estimated directly by

\[
\widehat f_n(\alpha)
=
\frac{\log N_n(\alpha,\Delta\alpha)}
{n\alpha}.
\]

For the minimal bRIFS with

\[
r\sim\mathrm{Uniform}(0,1),
\]

the analytic prediction is

\[
f(\alpha)
=
\frac{\log(m\alpha)+1-\alpha}{\alpha}.
\]

Example:

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

This produces:

- `Figure_1.pdf` — finite-depth contraction exponents and the depth-25 multifractal-spectrum estimate compared with the analytic prediction.
- `Figure_S1_Convergence.pdf` — finite-depth estimates at recursive depths 15, 20, and 25 compared with the analytic spectrum.

The same stochastic realization is used at every displayed depth.

---

### `spectrum_from_S.m`

Legacy scale-matrix implementation retained for reference.

This function analyzes the reconstructed scale matrix `S` using the earlier partition-sum and Legendre-transform workflow.

It is **not** the method used for the final numerical multifractal-spectrum figure in the manuscript.

The manuscript analysis should be reproduced with `spectrum_from_Leaves.m`.

## Reproducing the manuscript analysis from a new realization

From the repository root:

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

Because the bRIFS is stochastic, a newly generated realization will not exactly reproduce the saved manuscript realization. The saved data in the repository's `data/` directory should be used when exact reproduction of the manuscript numerical figure is desired.

## MATLAB requirements

The code was developed in MATLAB.

`BuildMultifractalTreeFn.m` currently uses:

- `parfor`
- `gpuArray`

and may therefore require the Parallel Computing Toolbox and a supported GPU for full simulation generation.

The saved manuscript data can be analyzed with `spectrum_from_Leaves.m` without rerunning the full bRIFS simulation.

## Important implementation note

Do not deduplicate the entries in `Leaves`.

Multiple descendant lineages may inherit the same contraction scale within a recursive replacement event. These are still distinct recursive lineages and must be counted separately in

\[
N_n(\alpha,\Delta\alpha).
\]

Removing duplicate scale trajectories changes lineage multiplicities and biases the finite-depth multifractal estimate.
