function [alphaCenters, fhatFinal, conv] = spectrum_from_Leaves(Leaves, m, mainOutFile, convOutFile, nBins, convDepths)
% SPECTRUM_FROM_LEAVES
% Direct finite-depth multifractal analysis of the bRIFS using the exact
% recursive lineage scales stored in Leaves.
%
% Each Leaves{n} must contain ONE scale value B_v for every recursive
% lineage/node v at recursive depth n, with multiplicity preserved.
%
% This matches the manuscript definitions:
%
%   alpha_n(v) = -log(B_v)/n
%
%   N_n(alpha,Delta alpha)
%       = number of depth-n lineages in an alpha bin
%
%   fhat_n(alpha)
%       = log N_n(alpha,Delta alpha)/(n*alpha)
%
% For the minimal bRIFS with r ~ Uniform(0,1), the analytic spectrum is
%
%   f(alpha) = [log(m*alpha) + 1 - alpha]/alpha.
%
% INPUTS
%   Leaves      : cell array; Leaves{n} contains all lineage scales at
%                 recursive depth n, including repeated equal scales for
%                 distinct sibling lineages.
%
%   m           : E[K], expected number of terminal descendants produced
%                 by one recursive replacement. For the manuscript model
%                 MaxOffspring=3, MaxGens=2, m ~= 1.657.
%
%   mainOutFile : optional filename for main two-panel figure.
%
%   convOutFile : optional filename for convergence figure.
%
%   nBins       : optional number of common alpha bins. Default = 30.
%
%   convDepths  : optional recursive depths for convergence comparison.
%                 Default = [15 20 numel(Leaves)], restricted to available
%                 depths and duplicates removed.
%
% OUTPUTS
%   alphaCenters : common alpha-bin centers
%   fhatFinal    : direct spectrum estimate at deepest available depth
%   conv         : struct with convergence depths and estimated spectra
%
% IMPORTANT
%   Do not deduplicate Leaves. Distinct sibling lineages may share exactly
%   the same inherited scale and must remain distinct in N_n(alpha).

    if nargin < 2 || isempty(m)
        error('Supply m = E[K]. For the manuscript model use approximately 1.657.');
    end
    if nargin < 3, mainOutFile = ''; end
    if nargin < 4, convOutFile = ''; end
    if nargin < 5 || isempty(nBins), nBins = 30; end
    if nargin < 6 || isempty(convDepths)
        convDepths = [15 20 numel(Leaves)];
    end

    if ~iscell(Leaves) || isempty(Leaves)
        error('Leaves must be a nonempty cell array.');
    end
    if ~isscalar(m) || ~isfinite(m) || m <= 1
        error('m must be a finite scalar greater than 1.');
    end
    if ~isscalar(nBins) || ~isfinite(nBins) || nBins < 5
        error('nBins must be a finite scalar >= 5.');
    end
    nBins = round(nBins);

    Ndepth = numel(Leaves);

    % ---------- Clean scale vectors and compute alpha_n ----------
    scales = cell(Ndepth,1);
    alphaByDepth = cell(Ndepth,1);
    lineageCounts = zeros(Ndepth,1);

    for n = 1:Ndepth
        b = double(Leaves{n}(:));
        b = b(isfinite(b) & b > 0);

        scales{n} = b;
        lineageCounts(n) = numel(b);

        if isempty(b)
            alphaByDepth{n} = [];
        else
            alphaByDepth{n} = -log(b) / n;
        end
    end

    validDepths = find(lineageCounts > 0);
    if isempty(validDepths)
        error('Leaves contains no positive finite scales.');
    end

    nFinal = validDepths(end);
    alphaFinal = alphaByDepth{nFinal};
    alphaFinal = alphaFinal(isfinite(alphaFinal) & alphaFinal > 0);

    if numel(alphaFinal) < 2
        error('Need at least two lineages at the deepest valid depth.');
    end

    % ---------- Convergence depths ----------
    convDepths = unique(round(convDepths(:).'), 'stable');
    convDepths = convDepths(convDepths >= 1 & convDepths <= Ndepth);
    convDepths = convDepths(lineageCounts(convDepths) > 0);

    if ~ismember(nFinal, convDepths)
        convDepths(end+1) = nFinal;
    end
    convDepths = unique(convDepths, 'stable');

    % ---------- Common alpha bins ----------
    % Use one common set of bins for all displayed depths. This avoids
    % changing Delta alpha with depth, which would confound convergence.
    pooled = [];
    for n = convDepths
        a = alphaByDepth{n};
        pooled = [pooled; a(isfinite(a) & a > 0)]; %#ok<AGROW>
    end

    alphaMin = min(pooled);
    alphaMax = max(pooled);

    if alphaMax <= alphaMin
        pad = max(1e-6, 0.05 * max(1, abs(alphaMin)));
        alphaMin = max(eps, alphaMin - pad);
        alphaMax = alphaMax + pad;
    end

    edges = linspace(alphaMin, alphaMax, nBins + 1);
    alphaCenters = (edges(1:end-1) + edges(2:end)) / 2;

    % ---------- Direct spectra at requested depths ----------
    fhatByDepth = nan(numel(convDepths), nBins);
    countsByDepth = zeros(numel(convDepths), nBins);

    for j = 1:numel(convDepths)
        n = convDepths(j);
        a = alphaByDepth{n};
        a = a(isfinite(a) & a > 0);

        counts = histcounts(a, edges);
        countsByDepth(j,:) = counts;

        good = counts > 0 & alphaCenters > 0;
        fhatByDepth(j,good) = ...
            log(counts(good)) ./ (n .* alphaCenters(good));
    end

    finalRow = find(convDepths == nFinal, 1, 'last');
    fhatFinal = fhatByDepth(finalRow,:);

    % ---------- Theory ----------
    theoryGridMin = max(eps, min(alphaMin, 0.05));
    theoryGridMax = max([alphaMax, 2, 2*log(m + 1)]);
    alphaTheory = linspace(theoryGridMin, theoryGridMax, 2500);

    theoryNumerator = log(m .* alphaTheory) + 1 - alphaTheory;
    fTheory = theoryNumerator ./ alphaTheory;
    attainable = theoryNumerator >= 0 & isfinite(fTheory);

    % ---------- Main figure ----------
    % Top: actual lineage-count distribution through recursive depth.
    % Bottom: final-depth direct spectrum + analytic prediction.
    allAlpha = vertcat(alphaByDepth{validDepths});
    allAlpha = allAlpha(isfinite(allAlpha) & allAlpha > 0);

    heatMin = min(allAlpha);
    heatMax = max(allAlpha);

    if heatMax <= heatMin
        pad = max(1e-6, 0.05 * max(1, abs(heatMin)));
        heatMin = max(eps, heatMin - pad);
        heatMax = heatMax + pad;
    end

    heatEdges = linspace(heatMin, heatMax, nBins + 1);
    heatCenters = (heatEdges(1:end-1) + heatEdges(2:end)) / 2;
    H = zeros(numel(validDepths), nBins);

    for j = 1:numel(validDepths)
        n = validDepths(j);
        a = alphaByDepth{n};
        a = a(isfinite(a) & a > 0);
        H(j,:) = histcounts(a, heatEdges);
    end

    f1 = figure('Color','w','Position',[60 60 950 700]);
    tl = tiledlayout(f1, 2, 1, ...
        'TileSpacing','compact', 'Padding','compact');

    ax1 = nexttile(tl,1);
    imagesc(ax1, heatCenters, validDepths, log1p(H));
    set(ax1,'YDir','normal');
    colormap(ax1, parula);
    cb = colorbar(ax1);
    ylabel(cb,'ln(count + 1)');
    xlabel(ax1,'\alpha_n = -ln(B_v)/n');
    ylabel(ax1,'recursive depth n');
    title(ax1,'Finite-depth contraction exponents across recursive lineages');
    box(ax1,'on');

    ax2 = nexttile(tl,2);
    hold(ax2,'on');

    goodFinal = isfinite(fhatFinal);
    plot(ax2, alphaCenters(goodFinal), fhatFinal(goodFinal), ...
        'o-', 'LineWidth',1.2, ...
        'DisplayName',sprintf('finite depth n = %d',nFinal));

    plot(ax2, alphaTheory(attainable), fTheory(attainable), ...
        '-', 'LineWidth',1.8, ...
        'DisplayName',sprintf('analytic uniform case, m = %.3f',m));

    grid(ax2,'on');
    box(ax2,'on');
    xlabel(ax2,'\alpha');
    ylabel(ax2,'f(\alpha)');
    title(ax2,'Direct finite-depth estimate of the multifractal spectrum');
    legend(ax2,'Location','best');

    sgtitle(tl,sprintf('bRIFS scale spectrum (m = %.3f)',m));

    if ~isempty(mainOutFile)
        saveFigure(f1, mainOutFile);
    end

    % ---------- Appendix convergence figure ----------
    f2 = figure('Color','w','Position',[100 100 900 560]);
    ax = axes(f2);
    hold(ax,'on');

    for j = 1:numel(convDepths)
        good = isfinite(fhatByDepth(j,:));
        plot(ax, alphaCenters(good), fhatByDepth(j,good), ...
            'o-', 'LineWidth',1.1, ...
            'DisplayName',sprintf('n = %d',convDepths(j)));
    end

    plot(ax, alphaTheory(attainable), fTheory(attainable), ...
        '-', 'LineWidth',2.0, ...
        'DisplayName','analytic spectrum');

    grid(ax,'on');
    box(ax,'on');
    xlabel(ax,'\alpha');
    ylabel(ax,'\hat{f}_n(\alpha)');
    title(ax,'Finite-depth convergence toward the analytic spectrum');
    legend(ax,'Location','best');

    if ~isempty(convOutFile)
        saveFigure(f2, convOutFile);
    end

    % ---------- Return convergence data ----------
    conv = struct();
    conv.depths = convDepths;
    conv.alphaCenters = alphaCenters;
    conv.counts = countsByDepth;
    conv.fhat = fhatByDepth;
    conv.alphaTheory = alphaTheory(attainable);
    conv.fTheory = fTheory(attainable);
    conv.lineageCounts = lineageCounts;

    % ---------- helper ----------
    function saveFigure(figHandle, filename)
        [~,~,ext] = fileparts(filename);

        switch lower(ext)
            case '.pdf'
                exportgraphics(figHandle, filename, 'ContentType','vector');
            case '.png'
                exportgraphics(figHandle, filename, 'Resolution',300);
            otherwise
                warning('Unknown extension; saving PNG.');
                exportgraphics(figHandle, [filename '.png'], 'Resolution',300);
        end
    end
end
