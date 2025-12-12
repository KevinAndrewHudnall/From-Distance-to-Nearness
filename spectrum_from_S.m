function [qVals, tau, alpha, falpha] = spectrum_from_S(S, qVals, outFile)
% SPECTRUM_FROM_S
%   Top   : ln-mass heatmap by depth with top-k overlay (visualization only)
%   Bottom: multifractal spectrum f(α) computed from RAW S (no renorm)
%
% INPUTS
%   S       : [Ndepth x Npaths] matrix, S(n,j) = s_j^(n) > 0 (scales)
%   qVals   : (optional) vector of q values, default linspace(-4,4,161)
%   outFile : (optional) filename to save figure ('.pdf' or '.png')
%
% OUTPUTS
%   qVals, tau(q), alpha(q), f(alpha)

    if nargin < 2 || isempty(qVals), qVals = linspace(-4,4,161); end
    if nargin < 3, outFile = ''; end

    % --- Clean S ---
    S = double(S);
    S(~isfinite(S) | S <= 0) = NaN;

    N = size(S,1);
    haveRow = sum(isfinite(S),2) >= 2;
    if nnz(haveRow) < 3
        error('Need ≥3 depths with ≥2 positive entries for a reliable spectrum.');
    end
    depths = find(haveRow).';
    D = numel(depths);

    % ---------- Build row-normalized masses (for HEATMAP only) ----------
    masses_by_depth = cell(D,1);
    for t = 1:D
        n = depths(t);
        row = S(n, :);
        row = row(isfinite(row));
        if isempty(row), row = NaN; end
        masses_by_depth{t} = row / nansum(row);
    end

    % ---------- Figure & layout ----------
    f = figure('Color','w','Position',[60 60 900 640]);
    tl = tiledlayout(f, 2, 1, 'TileSpacing','compact', 'Padding','compact');

    % ---------- Panel A: ln-mass heatmap + top-k overlay ----------
    axA = nexttile(tl, 1);  hold(axA,'on');
    W = 1024;                        % horizontal resolution of heatmap
    R = nan(D, W);
    eps0 = 1e-15;

    for t = 1:D
        p = masses_by_depth{t};
        R(t, :) = rebin_to_W(p, W);        % piecewise-constant rasterization
    end

    Rln = log(R + eps0);                    % natural log
    imagesc(axA, linspace(0,1,W), 0:D-1, Rln);
    set(axA,'YDir','normal');
    colormap(axA, parula);
    c = colorbar(axA); ylabel(c, 'ln(mass)');
    xlabel(axA, 'normalized position'); ylabel(axA, 'depth');
    title(axA, 'Mass across position & depth (natural log scale)');
    xlim(axA,[0 1]); ylim(axA,[0 D-1]);

    % overlay: top-k interval centers per depth
    k = 0;
    for t = 1:D
        p = masses_by_depth{t};
        M = numel(p); if M==0 || all(~isfinite(p)), continue; end
        [~, I] = maxk(p, min(k, M));
        cdfp = [0, cumsum(p(:).')];
        centers = (cdfp(I) + p(I)/2);
        y = (t-1) * ones(size(centers));
        plot(axA, centers, y, 'w.', 'MarkerSize', 10, 'LineWidth', 1.5);
    end

    % ---------- Panel B: MULTIFRACTAL SPECTRUM from RAW S ----------
    % representative scale ε_n: geometric mean of RAW S row
    eps_n = nan(D,1);
    for t = 1:D
        n = depths(t);
        rowS = S(n,:); rowS = rowS(isfinite(rowS));
        eps_n(t) = exp(mean(log(rowS)));
    end
    log_eps = log(eps_n);                     % natural log

    qVals = qVals(:).';
    Z = zeros(numel(qVals), D);
    for kq = 1:numel(qVals)
        q = qVals(kq);
        for t = 1:D
            n = depths(t);
            rowS = S(n,:); rowS = rowS(isfinite(rowS));
            Z(kq,t) = nansum(rowS.^q);        % <-- RAW S (no renorm)
        end
    end

    % Fit tau(q): slope of ln Z_n(q) vs ln ε_n
    tau = zeros(size(qVals));
    X = [log_eps(:), ones(numel(log_eps),1)];
    for kq = 1:numel(qVals)
        y = log(Z(kq,:)+eps).';
        pcoef = X \ y;                        % least-squares
        tau(kq) = pcoef(1);
    end

    % Legendre transform
    dq = gradient(qVals);
    alpha  = gradient(tau) ./ dq;
    falpha = qVals .* alpha - tau;

    axC = nexttile(tl, 2);
    plot(axC, alpha, falpha, 'o-', 'LineWidth', 1);
    grid(axC,'on'); box(axC,'on');
    xlabel(axC, '\alpha'); ylabel(axC, 'f(\alpha)');
    title(axC, 'Estimated multifractal spectrum f(\alpha)');

    % ---------- Save ----------
    if ~isempty(outFile)
        [~,~,ext] = fileparts(outFile);
        switch lower(ext)
            case '.pdf'
                set(f,'PaperPositionMode','auto'); print(f, outFile, '-dpdf','-r300');
            case '.png'
                exportgraphics(f, outFile, 'Resolution', 300);
            otherwise
                warning('Unknown extension; saving PNG.');
                exportgraphics(f, [outFile '.png'], 'Resolution', 300);
        end
    end

    % ---------- helper ----------
    function rowW = rebin_to_W(p, W)
        % piecewise-constant rebinning of masses p to length-W raster
        if isempty(p) || all(~isfinite(p))
            rowW = nan(1,W); return;
        end
        p = p(:).';
        p = p / sum(p);
        c = [0, cumsum(p)];
        x = (0.5:W-0.5)/W;
        [~, bins] = histc(x, c);
        bins(bins < 1) = 1; bins(bins > numel(p)) = numel(p);
        rowW = p(bins);
    end
end
