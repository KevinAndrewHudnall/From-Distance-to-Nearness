function [alphaCenters, fhat, alphaByDepth, mUsed] = spectrum_from_S(S, m, outFile, nBins, depthVals)
% SPECTRUM_FROM_S
%   Direct finite-depth estimate of the bRIFS multifractal spectrum from
%   inherited lineage scales S, matching the branching-random-walk theorem.
%
%   Top panel:
%       Distribution of finite-depth contraction exponents
%
%           alpha_n(v) = -log(B_v) / n
%
%       across recursive depth.
%
%   Bottom panel:
%       Direct finite-depth spectrum estimate
%
%           fhat_n(alpha) = log N_n(alpha,Delta alpha) / (n * alpha),
%
%       evaluated at the deepest available recursive depth, together with
%       the analytic spectrum for the minimal bRIFS with r ~ Uniform(0,1):
%
%           f(alpha) = [log(m*alpha) + 1 - alpha] / alpha.
%
% INPUTS
%   S         : [Ndepth x Npaths] matrix of inherited scales.
%               S(row,j) = B_v for lineage/node j at the corresponding
%               recursive depth. Missing/nonexistent entries may be NaN,
%               zero, or nonfinite.
%
%   m         : (optional) expected number of terminal descendants per
%               recursive replacement, m = E[K].
%               If empty, m is estimated from exponential growth in the
%               number of finite entries per recursive depth.
%
%   outFile   : (optional) filename for the figure ('.pdf' or '.png').
%
%   nBins     : (optional) number of alpha bins. Default = 30.
%
%   depthVals : (optional) recursive depth represented by each row of S.
%               Default = (1:size(S,1))'.
%               If row 1 is the root, pass (0:size(S,1)-1)'.
%               Rows with depth <= 0 are excluded from alpha calculations.
%
% OUTPUTS
%   alphaCenters : centers of alpha bins used for the final-depth estimate
%   fhat         : direct finite-depth estimate of f(alpha)
%   alphaByDepth : matrix, same size as S, containing -log(S)/depth
%   mUsed        : supplied m, or estimated m if m was empty
%
% NOTES
%   1. This function intentionally does NOT normalize S into probability
%      masses, compute partition sums Z_n(q), fit tau(q), or apply a
%      Legendre transform.
%
%   2. The direct estimator is a finite-depth approximation. Bins with very
%      few lineages are intrinsically noisy, especially near the boundaries
%      of the attainable alpha range.
%
%   3. The analytic overlay is specific to the minimal bRIFS with
%      r ~ Uniform(0,1). The empirical estimator itself is distribution-free.

    if nargin < 2
        m = [];
    end
    if nargin < 3 || isempty(outFile)
        outFile = '';
    end
    if nargin < 4 || isempty(nBins)
        nBins = 30;
    end
    if nargin < 5 || isempty(depthVals)
        depthVals = (1:size(S,1)).';
    end

    % ---------- Validate and clean inputs ----------
    S = double(S);
    S(~isfinite(S) | S <= 0) = NaN;

    depthVals = double(depthVals(:));
    if numel(depthVals) ~= size(S,1)
        error('depthVals must contain one recursive depth for each row of S.');
    end
    if any(~isfinite(depthVals))
        error('depthVals must contain only finite values.');
    end
    if ~isscalar(nBins) || ~isfinite(nBins) || nBins < 5
        error('nBins must be a finite scalar >= 5.');
    end
    nBins = round(nBins);

    lineageCounts = sum(isfinite(S), 2);

    % Need positive recursive depth and at least one valid lineage.
    validRows = depthVals > 0 & lineageCounts >= 1;
    if nnz(validRows) < 2
        error('Need at least two positive recursive depths containing valid scales.');
    end

    validIdx = find(validRows);
    validDepths = depthVals(validIdx);
    D = numel(validIdx);

    % ---------- Finite-depth contraction exponents ----------
    alphaByDepth = nan(size(S));

    for t = 1:D
        rowIdx = validIdx(t);
        n = depthVals(rowIdx);

        good = isfinite(S(rowIdx,:));
        alphaByDepth(rowIdx,good) = -log(S(rowIdx,good)) / n;
    end

    % ---------- Determine m ----------
    if isempty(m)
        % Estimate log(m) from growth in lineage number:
        % log L_n approximately n log(m) + constant.
        %
        % Use the deeper half of available rows to reduce early-depth
        % transients in a single realization.
        positiveCountRows = validRows & lineageCounts > 0;
        idxGrowth = find(positiveCountRows);

        if numel(idxGrowth) < 2
            error('Cannot estimate m: insufficient depths with positive lineage counts.');
        end

        firstFit = max(1, floor(numel(idxGrowth)/2));
        idxFit = idxGrowth(firstFit:end);

        x = depthVals(idxFit);
        y = log(lineageCounts(idxFit));

        if numel(unique(x)) < 2
            error('Cannot estimate m: depthVals do not contain enough distinct depths.');
        end

        p = polyfit(x, y, 1);
        mUsed = exp(p(1));

        if ~(isfinite(mUsed) && mUsed > 1)
            warning(['Estimated m = %.4g is not > 1. ' ...
                     'The uniform-case analytic spectrum will not be plotted.'], mUsed);
        end
        mWasEstimated = true;
    else
        if ~isscalar(m) || ~isfinite(m) || m <= 1
            error('m must be a finite scalar greater than 1.');
        end
        mUsed = double(m);
        mWasEstimated = false;
    end

    % ---------- Panel A data: alpha distribution through depth ----------
    allAlpha = alphaByDepth(validRows,:);
    allAlpha = allAlpha(isfinite(allAlpha) & allAlpha > 0);

    if isempty(allAlpha)
        error('No positive finite contraction exponents could be computed.');
    end

    alphaMinAll = min(allAlpha);
    alphaMaxAll = max(allAlpha);

    if alphaMaxAll <= alphaMinAll
        pad = max(1e-6, 0.05 * max(1, abs(alphaMinAll)));
        alphaMinAll = max(eps, alphaMinAll - pad);
        alphaMaxAll = alphaMaxAll + pad;
    end

    edgesAll = linspace(alphaMinAll, alphaMaxAll, nBins + 1);
    centersAll = (edgesAll(1:end-1) + edgesAll(2:end)) / 2;

    H = zeros(D, nBins);
    for t = 1:D
        rowIdx = validIdx(t);
        a = alphaByDepth(rowIdx,:);
        a = a(isfinite(a) & a > 0);
        H(t,:) = histcounts(a, edgesAll);
    end

    % ---------- Final-depth direct spectrum estimate ----------
    deepestRow = validIdx(end);
    nFinal = depthVals(deepestRow);

    alphaFinal = alphaByDepth(deepestRow,:);
    alphaFinal = alphaFinal(isfinite(alphaFinal) & alphaFinal > 0);

    if numel(alphaFinal) < 2
        error('Need at least two lineages at the deepest valid depth.');
    end

    alphaMin = min(alphaFinal);
    alphaMax = max(alphaFinal);

    if alphaMax <= alphaMin
        pad = max(1e-6, 0.05 * max(1, abs(alphaMin)));
        alphaMin = max(eps, alphaMin - pad);
        alphaMax = alphaMax + pad;
    end

    edges = linspace(alphaMin, alphaMax, nBins + 1);
    alphaCenters = (edges(1:end-1) + edges(2:end)) / 2;

    counts = histcounts(alphaFinal, edges);

    fhat = nan(size(alphaCenters));
    nonempty = counts > 0 & alphaCenters > 0;
    fhat(nonempty) = log(counts(nonempty)) ./ (nFinal .* alphaCenters(nonempty));

    % ---------- Figure ----------
    f = figure('Color','w','Position',[60 60 900 680]);
    tl = tiledlayout(f, 2, 1, 'TileSpacing','compact', 'Padding','compact');

    % Panel A: log lineage counts across alpha and recursive depth
    axA = nexttile(tl, 1);
    imagesc(axA, centersAll, validDepths, log1p(H));
    set(axA, 'YDir', 'normal');
    colormap(axA, parula);
    c = colorbar(axA);
    ylabel(c, 'ln(count + 1)');
    xlabel(axA, '\alpha_n = -ln(B_v)/n');
    ylabel(axA, 'recursive depth n');
    title(axA, 'Finite-depth contraction exponents across recursive lineages');
    box(axA, 'on');

    % Panel B: direct finite-depth estimate + analytic uniform spectrum
    axB = nexttile(tl, 2);
    hold(axB, 'on');

    plotMask = isfinite(fhat);
    plot(axB, alphaCenters(plotMask), fhat(plotMask), 'o-', ...
        'LineWidth', 1.2, 'DisplayName', sprintf('finite depth n = %g', nFinal));

    if isfinite(mUsed) && mUsed > 1
        % Plot analytic spectrum over a broad positive range, then retain
        % only the attainable domain where log(m*alpha)+1-alpha >= 0.
        gridMin = max(eps, min(alphaMin, 0.05));
        gridMax = max([alphaMax, 2, 2*log(mUsed + 1)]);
        alphaTheory = linspace(gridMin, gridMax, 2000);

        numerator = log(mUsed .* alphaTheory) + 1 - alphaTheory;
        fTheory = numerator ./ alphaTheory;

        attainable = numerator >= 0 & isfinite(fTheory);
        plot(axB, alphaTheory(attainable), fTheory(attainable), '-', ...
            'LineWidth', 1.8, ...
            'DisplayName', sprintf('analytic uniform case, m = %.4g', mUsed));
    end

    grid(axB, 'on');
    box(axB, 'on');
    xlabel(axB, '\alpha');
    ylabel(axB, 'f(\alpha)');
    title(axB, 'Direct finite-depth estimate of the multifractal spectrum');
    legend(axB, 'Location', 'best');

    if mWasEstimated
        sgtitle(tl, sprintf('bRIFS scale spectrum (m estimated as %.4g)', mUsed));
    else
        sgtitle(tl, sprintf('bRIFS scale spectrum (m = %.4g)', mUsed));
    end

    % ---------- Save ----------
    if ~isempty(outFile)
        [~,~,ext] = fileparts(outFile);

        switch lower(ext)
            case '.pdf'
                exportgraphics(f, outFile, 'ContentType', 'vector');
            case '.png'
                exportgraphics(f, outFile, 'Resolution', 300);
            otherwise
                warning('Unknown extension; saving PNG.');
                exportgraphics(f, [outFile '.png'], 'Resolution', 300);
        end
    end
end
