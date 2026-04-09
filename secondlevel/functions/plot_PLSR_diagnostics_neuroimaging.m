function ROI_table = plot_PLSR_diagnostics_neuroimaging(results, XL, roiNames, atlasFile, varargin)
%
% Plot/diagnose PLSR neuroimaging results (ROI/parcellation).
%
% This function provides generic visualization and diagnostics for PLSR
% results produced by pipelines such as PLSR_neuroimaging_pipeline.
% It creates:
%   - ROI table export (VIP + stabilityZ + robust flag, plus optional beta/stability fields)
%   - VIP vs stabilityZ scatter with robust ROI labeling
%   - NIfTI maps: VIP, stabilityZ, LV loadings (raw + thresholded)
%   - Multi-slice visualization (LV / VIP / stabilityZ) with ROI labels
%   - Optional: top-K mean weights table/bar plot (if results.meanFeatureWeight or results.meanBeta exists)
%   - Optional: feature stability stem plot (if results.featureStability exists)
%   - Optional: predicted vs observed scatter (if held-out CV predictions are stored)
%   - Optional: permutation Q2 histogram with observed Q2 overlay
%   - Optional: bootstrap Q2 histogram with observed Q2 overlay
%   - Optional: learning curve plot
%
% USAGE
%   ROI_table = plot_PLSR_diagnostics_neuroimaging(results, XL, roiNames, atlasFile)
%   ROI_table = plot_PLSR_diagnostics_neuroimaging(results, XL, roiNames, atlasFile, 'Name',Value,...)
%
% INPUTS
%   results struct
%        Must contain:
%          results.VIP        [p x 1] VIP scores from the final PLS model
%          results.stabilityZ [p x 1] stability statistic (meanBeta ./ sdBeta)
%        Optional fields used if present:
%          results.finalXLoadings     [p x finalLV] used if XL is empty
%          results.finalLV            scalar       clamps LV selection if present
%          results.meanFeatureWeight  [p x 1] for top-K weight table/bar plot
%          results.meanBeta           [p x 1] preferred over meanFeatureWeight if present
%          results.selectionFrequency [p x 1] for ROI table / top-K table
%          results.signStability      [p x 1] for ROI table / top-K table
%          results.featureStability   [p x 1] for stability stem plot
%          results.cvObserved         [N x 1] observed Y values from outer CV (optional)
%          results.cvPredicted        [N x 1] predicted Y values from outer CV (optional)
%          results.cvRepeatID         [N x 1] repeat index for each held-out prediction (optional)
%          results.cvSubjectID        [N x 1] subject index for each held-out prediction (optional)
%          results.allpermQ2          [nPerm x 1] permutation distribution (optional)
%          results.Q2                 scalar observed nested CV Q2 (optional but recommended)
%          results.allbootQ2          [nBoot x 1] bootstrap OOB Q2 distribution (optional)
%          results.Q2_CI              [1 x 2] bootstrap CI (optional)
%          results.learningSizes      vector sample sizes (optional)
%          results.learningQ2         vector Q2 values (optional)
%          results.RMSE               scalar mean RMSE from nested CV (optional)
%
%   XL   [p x L] numeric
%        X-loadings from plsregress (columns = LVs). Used to create the LV
%        brain map. If empty, the function tries results.finalXLoadings.
%
%   roiNames  cellstr/string (length p)
%        ROI/feature labels. If empty, defaults to 'ROI_001'...'ROI_p'.
%
%   atlasFile char/string
%        Path to labeled atlas NIfTI. Voxels should be labeled with integers
%        1..p matching the feature order in X/results/XL (0 = background).
%        If empty, atlas-based NIfTI mapping and slice visualization are skipped.
%
% NAME–VALUE OPTIONS
%   'TopN'        (20)   Number of ROI-table rows printed to console (CSV exports all rows).
%   'LV'          (1)    LV index to visualize (clamped to valid range).
%   'OutPrefix'   ('PLSR')  Prefix for exported CSV/NIfTI files.
%   'VIP_thresh'  (1)    VIP threshold for robust contributor definition/labeling.
%   'stab_thresh' (2)    stabilityZ threshold for robust contributor definition/labeling.
%   'MapPrctile'  (70)   Threshold LV map by abs percentile (keeps top 100-MapPrctile %).
%   'RelaxIfEmpty'(true) If no ROI passes thresholds, relax thresholds (75th percentiles)
%                        for visualization/labeling only.
%   'UnderlayFile' ('')   Optional structural underlay NIfTI for the multi-slice figure.
%                         Must be in the same voxel space/dimensions as atlasFile.
%
% OUTPUT
%   ROI_table table
%        Core columns:
%          ROI, VIP, stabilityZ, RobustContributor
%        Additional columns included if available:
%          meanBeta, signStability, selectionFrequency, featureStability
%        Also written to: <OutPrefix>_ROI_VIP_stability.csv
%
% FILES WRITTEN (requires SPM for NIfTI I/O)
%   <OutPrefix>_ROI_VIP_stability.csv
%   <OutPrefix>_VIP_map.nii
%   <OutPrefix>_stabilityZ_map.nii
%   <OutPrefix>_LV<LV>_map.nii
%   <OutPrefix>_LV<LV>_map_thresh.nii
%
% FIGURES CREATED (interactive)
%   1) VIP vs stabilityZ scatter (robust ROIs labeled)
%   2) Multi-slice 3×N panel: LV map / VIP map / stabilityZ map (robust ROIs labeled)
%   3) (optional) Top-K mean beta/weight bar plot
%   4) (optional) Feature stability stem plot
%   5) (optional) Predicted vs observed scatter
%      - held-out predictions stacked across repeated outer CV
%      - identity line
%      - OLS fitted line
%      - robust regression line
%      - optional density underlay
%      - optional coloring by subject ID across repeats
%      - optional subject trajectories across repeats
%   6) (optional) Permutation Q2 histogram
%   7) (optional) Bootstrap Q2 histogram
%   8) (optional) Learning curve
%
% NOTES
%   - Robust contributors are defined as VIP>VIP_thresh AND |stabilityZ|>stab_thresh.
%   - Correct atlas/feature alignment is essential: atlas labels 1..p must match
%     feature order in X/results/XL.
%   - For the predicted-vs-observed plot, it is best if results.cvObserved and
%     results.cvPredicted contain outer-fold held-out values only.
%   - If results.cvSubjectID is available, the predicted-vs-observed plot can
%     color repeated predictions by subject across CV repetitions, which helps
%     visualize subject-specific prediction bias or instability.
%   - The predicted-vs-observed figure is intended primarily as a diagnostic
%     visualization; quantitative inference should rely on results.Q2,
%     permutation testing, bootstrap Q2, and related fold-wise metrics from the pipeline.
%
% DEPENDENCIES
%   Requires SPM (spm_vol, spm_read_vols, spm_write_vol) for NIfTI mapping.
%   The robust fitted line in the predicted-vs-observed panel uses robustfit
%   if available.
%
% See also: plsregress, robustfit, spm_vol, spm_read_vols, spm_write_vol

% --------------------------
% Options
% --------------------------
p_in = inputParser;
addParameter(p_in,'TopN',20,@(x) isnumeric(x) && isscalar(x));
addParameter(p_in,'LV',1,@(x) isnumeric(x) && isscalar(x));
addParameter(p_in,'OutPrefix','PLSR',@(x) ischar(x) || isstring(x));
addParameter(p_in,'VIP_thresh',1,@(x) isnumeric(x) && isscalar(x));
addParameter(p_in,'stab_thresh',2,@(x) isnumeric(x) && isscalar(x));
addParameter(p_in,'MapPrctile',70,@(x) isnumeric(x) && isscalar(x));
addParameter(p_in,'RelaxIfEmpty',true,@(x) islogical(x) && isscalar(x));
addParameter(p_in,'UnderlayFile','',@(x) ischar(x) || isstring(x));
parse(p_in,varargin{:});

TopN         = p_in.Results.TopN;
LV           = p_in.Results.LV;
OutPrefix    = char(p_in.Results.OutPrefix);
VIP_thresh   = p_in.Results.VIP_thresh;
stab_thresh  = p_in.Results.stab_thresh;
MapPrctile   = p_in.Results.MapPrctile;
RelaxIfEmpty = p_in.Results.RelaxIfEmpty;
UnderlayFile = char(p_in.Results.UnderlayFile);

% --------------------------
% Checks / fallbacks
% --------------------------
if ~isfield(results,'VIP') || isempty(results.VIP)
    error('results.VIP missing');
end
if ~isfield(results,'stabilityZ') || isempty(results.stabilityZ)
    error('results.stabilityZ missing');
end

p = length(results.VIP);

if nargin < 2 || isempty(XL)
    if isfield(results,'finalXLoadings') && ~isempty(results.finalXLoadings)
        XL = results.finalXLoadings;
    else
        error('XL empty and results.finalXLoadings not found.');
    end
end

if nargin < 3 || isempty(roiNames)
    roiNames = arrayfun(@(i) sprintf('ROI_%03d',i), 1:p, 'UniformOutput', false);
end
roiNames = cellstr(roiNames(:));

if nargin < 4
    atlasFile = [];
end

% Robust LV selection/clamping
LV = max(1, min(LV, size(XL,2)));
% if isfield(results,'finalLV') && ~isempty(results.finalLV)
%     LV = max(1, min(LV, results.finalLV));
%     LV = max(1, min(LV, size(XL,2)));
% end

% --------------------------
% 1. ROI table
% --------------------------
isRobust = (results.VIP(:) > VIP_thresh) & (abs(results.stabilityZ(:)) > stab_thresh);

% Relax thresholds for visualization only if nothing passes
VIP_thresh_vis  = VIP_thresh;
stab_thresh_vis = stab_thresh;
if RelaxIfEmpty && ~any(isRobust)
    VIP_thresh_vis  = prctile(results.VIP(:),75);
    stab_thresh_vis = prctile(abs(results.stabilityZ(:)),75);
    isRobust = (results.VIP(:) > VIP_thresh_vis) & (abs(results.stabilityZ(:)) > stab_thresh_vis);
end

% Optional columns
if isfield(results,'meanBeta') && ~isempty(results.meanBeta)
    meanBeta = results.meanBeta(:);
elseif isfield(results,'meanFeatureWeight') && ~isempty(results.meanFeatureWeight)
    meanBeta = results.meanFeatureWeight(:);
else
    meanBeta = nan(p,1);
end

if isfield(results,'signStability') && ~isempty(results.signStability)
    signStability = results.signStability(:);
else
    signStability = nan(p,1);
end

ROI_table = table(roiNames, results.VIP(:), results.stabilityZ(:), meanBeta, ...
    signStability, isRobust(:), ...
    'VariableNames', {'ROI','VIP','stabilityZ','meanBeta','signStability', ...
    'RobustContributor'});

% Sort primarily by robustness-relevant quantities
[~, sortIdx] = sortrows([double(isRobust(:)) results.VIP(:) abs(results.stabilityZ(:)) abs(meanBeta)], ...
    [-1 -2 -3 -4]);
ROI_table = ROI_table(sortIdx,:);

disp(ROI_table(1:min(TopN,height(ROI_table)),:))
writetable(ROI_table, sprintf('%s_ROI_VIP_stability.csv',OutPrefix));

% --------------------------
% 2. Scatter VIP vs stabilityZ
% --------------------------
figure('Color','w','Units','pixels','Position',[100 100 1000 800]); hold on;
scatter(results.VIP(:), results.stabilityZ(:), 60, 'filled');
xlabel('VIP','FontSize',14);
ylabel('stabilityZ','FontSize',14);
title('ROI importance vs stability (PLSR)','FontSize',16);
grid on;
set(gca,'FontSize',14,'LineWidth',1.2);
xline(VIP_thresh_vis,'r--','LineWidth',1.2);
yline(stab_thresh_vis,'r--','LineWidth',1.2);
yline(-stab_thresh_vis,'r--','LineWidth',1.2);

robustIdx = find(isRobust);
for ii = 1:length(robustIdx)
    idx = robustIdx(ii);
    text(results.VIP(idx)+0.02, results.stabilityZ(idx), roiNames{idx}, ...
        'FontSize',11, 'Interpreter','none');
end

set(gcf,'Name','VIP_vs_stability','NumberTitle','off');

% --------------------------
% 3. Atlas load + sanity checks
% --------------------------
doAtlas = ~isempty(atlasFile);
if doAtlas
    V = spm_vol(atlasFile);
    atlasData = spm_read_vols(V);

    labels = unique(atlasData(:));
    labels(labels==0 | isnan(labels)) = [];
    if ~isempty(labels) && max(labels) < p
        warning('Atlas max label (%d) < number of features p (%d). Mapping may be wrong.', max(labels), p);
    end

    % Prefer slices that contain atlas content
    zHasData = squeeze(any(any(atlasData>0,1),2));
    zList = find(zHasData);
    if numel(zList) >= 5
        slice_idx = zList(round(linspace(1,numel(zList),5)));
    else
        slice_idx = round(linspace(1,size(atlasData,3),5));
    end
else
    atlasData = [];
    V = [];
    slice_idx = [];
end

% Optional underlay
hasUnderlay = doAtlas && ~isempty(UnderlayFile) && exist(UnderlayFile,'file');
if hasUnderlay
    Vu = spm_vol(UnderlayFile);
    underlayData = spm_read_vols(Vu);

    if ~isequal(size(underlayData), size(atlasData))
        error('Underlay and atlas must have the same dimensions.');
    end
else
    underlayData = [];
end

% --------------------------
% 4. VIP & stabilityZ maps
% --------------------------
if doAtlas
    VIPmap  = zeros(size(atlasData));
    stabMap = zeros(size(atlasData));

    for i = 1:p
        mask = atlasData == i;
        VIPmap(mask)  = results.VIP(i);
        stabMap(mask) = results.stabilityZ(i);
    end

    Vvip = V;   Vvip.fname  = sprintf('%s_VIP_map.nii',OutPrefix);
    Vstab = V;  Vstab.fname = sprintf('%s_stabilityZ_map.nii',OutPrefix);
    spm_write_vol(Vvip, VIPmap);
    spm_write_vol(Vstab, stabMap);
end

% --------------------------
% 5. LV map (raw + thresholded)
% --------------------------
LVweights = XL(:,LV);
den = max(abs(LVweights));
if den==0 || isnan(den) || isinf(den)
    den = 1;
end
LVweights = LVweights / den;

if doAtlas
    LVmap = zeros(size(atlasData));
    for i = 1:p
        mask = atlasData == i;
        LVmap(mask) = LVweights(i);
    end

    Vlv = V; Vlv.fname = sprintf('%s_LV%d_map.nii',OutPrefix,LV);
    spm_write_vol(Vlv, LVmap);

    th = prctile(abs(LVweights), MapPrctile);
    LVmap_thresh = zeros(size(atlasData));
    for i = 1:p
        mask = atlasData == i;
        val = LVweights(i);
        if abs(val) < th
            val = 0;
        end
        LVmap_thresh(mask) = val;
    end
    VlvT = V; VlvT.fname = sprintf('%s_LV%d_map_thresh.nii',OutPrefix,LV);
    spm_write_vol(VlvT, LVmap_thresh);
end

% --------------------------
% 6. Multi-slice figure: LV / VIP / stabilityZ
% --------------------------
if doAtlas
    % Recompute slice selection here so we can avoid extreme inferior/superior slices
    zHasData = squeeze(any(any(atlasData>0,1),2));
    zList = find(zHasData);

    if numel(zList) >= 5
        % Trim lower and upper extremes a bit
        lo = max(1, round(0.20 * numel(zList)));
        hi = min(numel(zList), round(0.80 * numel(zList)));
        keepIdx = zList(lo:hi);

        if numel(keepIdx) >= 5
            slice_idx = keepIdx(round(linspace(1, numel(keepIdx), 5)));
        else
            slice_idx = zList(round(linspace(1, numel(zList), 5)));
        end
    else
        slice_idx = round(linspace(1,size(atlasData,3),5));
    end

    VIPmap_vis = VIPmap;
    VIPmap_vis(VIPmap_vis < VIP_thresh_vis) = 0;

    stabMap_vis = stabMap;
    stabMap_vis(abs(stabMap_vis) < stab_thresh_vis) = 0;

    hFig = figure('Color','k','Units','pixels','Position',[50 50 2200 1200]);

    nSlices = length(slice_idx);
    baseAxes = gobjects(3,nSlices);

    for s = 1:nSlices
        z = slice_idx(s);

        % ---------- Row 1: LV ----------
        axBase = subplot(3,nSlices,s);
        baseAxes(1,s) = axBase;
        hold(axBase,'on');

        if hasUnderlay
            bg = orient_for_display(underlayData(:,:,z));
            bg = bg - min(bg(:));
            if max(bg(:)) > 0
                bg = bg ./ max(bg(:));
            end
        else
            bg = orient_for_display(double(atlasData(:,:,z) > 0));
        end

        imagesc(bg,'Parent',axBase);
        axis(axBase,'image');
        axis(axBase,'off');
        colormap(axBase,gray);

        pos = get(axBase,'Position');
        axOverlay = axes('Position',pos,'Color','none');
        ov = rot90(orient_for_display(LVmap(:,:,z)), 2);   % 180° flip for overlay only
        hOv = imagesc(ov,'Parent',axOverlay);
        axis(axOverlay,'image');
        axis(axOverlay,'off');
        colormap(axOverlay,jet);

        nz = ov(abs(ov) > 0);
        if ~isempty(nz)
            cmax = max(abs(nz));
            if cmax == 0 || isnan(cmax) || isinf(cmax), cmax = 1; end
            caxis(axOverlay,[-cmax cmax]);
        end

        set(hOv,'AlphaData', 0.75 * double(abs(ov) > 0));

        % ---------- Row 2: VIP ----------
        axBase = subplot(3,nSlices,s+nSlices);
        baseAxes(2,s) = axBase;
        hold(axBase,'on');

        if hasUnderlay
            bg = orient_for_display(underlayData(:,:,z));
            bg = bg - min(bg(:));
            if max(bg(:)) > 0
                bg = bg ./ max(bg(:));
            end
        else
            bg = orient_for_display(double(atlasData(:,:,z) > 0));
        end

        imagesc(bg,'Parent',axBase);
        axis(axBase,'image');
        axis(axBase,'off');
        colormap(axBase,gray);

        pos = get(axBase,'Position');
        axOverlay = axes('Position',pos,'Color','none');
        ov = rot90(orient_for_display(VIPmap_vis(:,:,z)), 2);   % 180° flip for overlay only
        hOv = imagesc(ov,'Parent',axOverlay);
        axis(axOverlay,'image');
        axis(axOverlay,'off');
        colormap(axOverlay,hot);

        nz = ov(ov > 0);
        if ~isempty(nz)
            cmin = min(nz);
            cmax = max(nz);
            if cmax <= cmin
                cmax = cmin + eps;
            end
            caxis(axOverlay,[cmin cmax]);
        end

        set(hOv,'AlphaData', 0.75 * double(ov > 0));

        % ---------- Row 3: stabilityZ ----------
        axBase = subplot(3,nSlices,s+2*nSlices);
        baseAxes(3,s) = axBase;
        hold(axBase,'on');

        if hasUnderlay
            bg = orient_for_display(underlayData(:,:,z));
            bg = bg - min(bg(:));
            if max(bg(:)) > 0
                bg = bg ./ max(bg(:));
            end
        else
            bg = orient_for_display(double(atlasData(:,:,z) > 0));
        end

        imagesc(bg,'Parent',axBase);
        axis(axBase,'image');
        axis(axBase,'off');
        colormap(axBase,gray);

        pos = get(axBase,'Position');
        axOverlay = axes('Position',pos,'Color','none');
        ov = rot90(orient_for_display(stabMap_vis(:,:,z)), 2);   % 180° flip for overlay only
        hOv = imagesc(ov,'Parent',axOverlay);
        axis(axOverlay,'image');
        axis(axOverlay,'off');
        colormap(axOverlay,jet);

        nz = ov(abs(ov) > 0);
        if ~isempty(nz)
            cmax = max(abs(nz));
            if cmax == 0 || isnan(cmax) || isinf(cmax), cmax = 1; end
            caxis(axOverlay,[-cmax cmax]);
        end

        set(hOv,'AlphaData', 0.75 * double(abs(ov) > 0));
    end

    % Add ROI labels onto the visible base axes
    for s = 1:nSlices
        z = slice_idx(s);

        for row = 1:3
            axBase = baseAxes(row,s);
            hold(axBase,'on');

            for ii = 1:length(robustIdx)
                roiID = robustIdx(ii);

                % IMPORTANT: label mask must follow overlay orientation
                mask = rot90(orient_for_display(atlasData(:,:,z) == roiID), 2);
                [yy,xx] = find(mask);

                % only label ROIs with enough visible voxels in this slice
                if numel(xx) >= 8
                    yPlot = size(mask,1) - median(yy) + 1;
                    text(axBase, median(xx), yPlot, roiNames{roiID}, ...
                        'Color','w', ...
                        'FontSize',8, ...
                        'FontWeight','bold', ...
                        'HorizontalAlignment','center', ...
                        'VerticalAlignment','middle', ...
                        'Interpreter','none', ...
                        'BackgroundColor','k', ...
                        'Margin',0.5);
                end
            end
        end
    end

    % Column headers: show real z in mm once on top row only
    for s = 1:nSlices
        ax = baseAxes(1,s);
        nx = size(atlasData,1);
        ny = size(atlasData,2);
        xyz_mm = V.mat * [nx/2; ny/2; slice_idx(s); 1];
        z_mm = xyz_mm(3);

        title(ax, sprintf('z = %.0f mm', z_mm), ...
            'Color','w', 'FontSize',16, 'FontWeight','bold');
    end

    % Row labels once at left
    rowNames = {'LV','VIP','stabZ'};
    for row = 1:3
        ax = baseAxes(row,1);
        text(ax, -0.10, 0.5, rowNames{row}, ...
            'Units','normalized', ...
            'Color','w', ...
            'FontSize',16, ...
            'FontWeight','bold', ...
            'HorizontalAlignment','right', ...
            'VerticalAlignment','middle', ...
            'Interpreter','none');
    end

    set(gcf,'Name','PLSR_multislice','NumberTitle','off');
end

% --------------------------
% 7. Top-N weights + stability (if present)
% --------------------------
if isfield(results,'meanFeatureWeight') && ~isempty(results.meanFeatureWeight)
    mf = results.meanFeatureWeight(:);
    [~,idx] = sort(abs(mf),'descend');

    if isfield(results,'signStability') && ~isempty(results.signStability)
        ss = results.signStability(:);
    else
        ss = nan(p,1);
    end

    weight_table = table(roiNames(idx(1:TopN)), mf(idx(1:TopN)), ss(idx(1:TopN)), ...
        'VariableNames', {'feature label','mean weight','sign stability'});
    disp(weight_table);
    writetable(weight_table, sprintf('%s_top%d_weights.csv',OutPrefix,TopN));

    figure('Color','w','Units','pixels','Position',[100 100 1400 800]);
    bar(mf(idx(1:TopN)));
    xticks(1:TopN);
    xticklabels(roiNames(idx(1:TopN)));
    set(gca,'TickLabelInterpreter','none');
    xtickangle(45);
    ylabel('Mean Weight','FontSize',14);
    title(sprintf('Top %d Features (mean weight)',TopN),'FontSize',16);
    set(gca,'FontSize',14,'LineWidth',1.2);
    set(gcf,'Name','TopN_weights','NumberTitle','off');
end

if isfield(results,'signStability') && ~isempty(results.signStability)
    figure('Color','w','Units','pixels','Position',[100 100 1400 700]);
    set(gcf,'Name','PLSR_sign_stability','NumberTitle','off');

    stem(results.signStability(:),'LineWidth',1.2);
    xlabel('Feature Index','FontSize',14);
    ylabel('Sign stability proportion','FontSize',14);
    title('Feature Sign Stability','FontSize',16);
    set(gca,'FontSize',14,'LineWidth',1.2);
end

% --------------------------
% 8. Predicted vs observed (if available)
% --------------------------
if isfield(results,'cvObserved') && isfield(results,'cvPredicted') && ...
        ~isempty(results.cvObserved) && ~isempty(results.cvPredicted)

    yObs = results.cvObserved(:);
    yPred = results.cvPredicted(:);

    hasRepeat = isfield(results,'cvRepeatID') && ~isempty(results.cvRepeatID);
    hasSubj   = isfield(results,'cvSubjectID') && ~isempty(results.cvSubjectID);

    if hasRepeat
        repID = results.cvRepeatID(:);
    else
        repID = nan(size(yObs));
    end

    if hasSubj
        subjID = results.cvSubjectID(:);
    else
        subjID = nan(size(yObs));
    end

    good = isfinite(yObs) & isfinite(yPred);
    if hasRepeat
        good = good & isfinite(repID);
    end
    if hasSubj
        good = good & isfinite(subjID);
    end

    yObs = yObs(good);
    yPred = yPred(good);
    if hasRepeat
        repID = repID(good);
    end
    if hasSubj
        subjID = subjID(good);
    end

    if ~isempty(yObs)

        figure('Color','w','Units','pixels','Position',[100 100 1000 800]); hold on;
        set(gcf,'Name','PLSR_predicted_vs_observed','NumberTitle','off');

        % ---------------------------------
        % Density underlay (2D histogram)
        % ---------------------------------
        showDensity = true;
        if showDensity && numel(yObs) >= 10

            nBins = max(20, round(sqrt(numel(yObs))));

            xMin = min(yObs);
            xMax = max(yObs);
            yMin = min(yPred);
            yMax = max(yPred);

            if xMin == xMax
                xMin = xMin - 0.5;
                xMax = xMax + 0.5;
            end
            if yMin == yMax
                yMin = yMin - 0.5;
                yMax = yMax + 0.5;
            end

            xEdges = linspace(xMin, xMax, nBins+1);
            yEdges = linspace(yMin, yMax, nBins+1);

            N = histcounts2(yObs, yPred, xEdges, yEdges);

            xCenters = (xEdges(1:end-1) + xEdges(2:end)) / 2;
            yCenters = (yEdges(1:end-1) + yEdges(2:end)) / 2;

            imagesc(xCenters, yCenters, N');
            set(gca, 'YDir', 'normal');

            hImg = findobj(gca, 'Type', 'image');
            if ~isempty(hImg)
                alphaData = double(N' > 0);
                alphaData(alphaData > 0) = 0.35;
                set(hImg(1), 'AlphaData', alphaData);
            end

            colormap(parula);
            colorbar;
        end

        % ---------------------------------
        % Scatter colored by subject ID if available
        % ---------------------------------
        if hasSubj

            uSubj = unique(subjID);
            nSubj = numel(uSubj);

            % Use a repeatable colormap that can handle many subjects
            cmap = lines(max(nSubj,7));
            if nSubj > size(cmap,1)
                cmap = hsv(nSubj);
            else
                cmap = cmap(1:nSubj,:);
            end

            for ii = 1:nSubj
                idx = subjID == uSubj(ii);

                scatter(yObs(idx), yPred(idx), 30, ...
                    'MarkerFaceColor', cmap(ii,:), ...
                    'MarkerEdgeColor', 'k', ...
                    'MarkerFaceAlpha', 0.75, ...
                    'MarkerEdgeAlpha', 0.4, ...
                    'HandleVisibility','off');
            end

            % Optional subject trajectories across repeats
            if hasRepeat
                for ii = 1:nSubj
                    idx = subjID == uSubj(ii);

                    if sum(idx) > 1
                        % sort within subject by repeat for cleaner trajectories
                        [~,ord] = sort(repID(idx));
                        xLine = yObs(idx);
                        yLine = yPred(idx);
                        xLine = xLine(ord);
                        yLine = yLine(ord);

                        plot(xLine, yLine, '-', ...
                            'Color', cmap(ii,:), ...
                            'LineWidth', 0.5, ...
                            'HandleVisibility','off');
                    end
                end
            end

            % Dummy legend handle
            scatter(nan, nan, 35, 'filled', ...
                'MarkerFaceColor', [0.4 0.4 0.4], ...
                'MarkerEdgeColor', 'k', ...
                'DisplayName', 'Subjects (colored by ID)');

        else
            scatter(yObs, yPred, 40, 'filled', ...
                'MarkerFaceAlpha', 0.8, ...
                'DisplayName', 'Subjects');
        end

        % ---------------------------------
        % Identity line
        % ---------------------------------
        mn = min([yObs; yPred]);
        mx = max([yObs; yPred]);
        if isfinite(mn) && isfinite(mx) && mn < mx
            plot([mn mx], [mn mx], 'k--', 'LineWidth', 1.5, ...
                'DisplayName', 'Identity line');
        end

        % ---------------------------------
        % OLS regression line
        % ---------------------------------
        bOLS = [];
        if numel(yObs) >= 2 && std(yObs) > 0
            bOLS = polyfit(yObs, yPred, 1);
            xfit = linspace(min(yObs), max(yObs), 100);
            yfitOLS = polyval(bOLS, xfit);
            plot(xfit, yfitOLS, 'r-', 'LineWidth', 2, ...
                'DisplayName', 'OLS fit');
        end

        % ---------------------------------
        % Robust regression line
        % ---------------------------------
        robustSlope = NaN;
        robustIntercept = NaN;
        try
            brob = robustfit(yObs, yPred);   % [intercept; slope]
            robustIntercept = brob(1);
            robustSlope = brob(2);

            xfitR = linspace(min(yObs), max(yObs), 100);
            yfitR = robustIntercept + robustSlope * xfitR;
            plot(xfitR, yfitR, 'm-', 'LineWidth', 2, ...
                'DisplayName', 'Robust fit');
        catch
            % robustfit unavailable or failed; skip silently
        end

        % ---------------------------------
        % Correlation
        % ---------------------------------
        rPlot = NaN;
        if numel(yObs) > 1 && std(yObs) > 0 && std(yPred) > 0
            C = corrcoef(yObs, yPred);
            rPlot = C(1,2);
        end

        % ---------------------------------
        % Display-only R^2
        % ---------------------------------
        R2plot = NaN;
        denomPlot = sum((yObs - mean(yObs)).^2);
        if denomPlot > 0
            R2plot = 1 - sum((yObs - yPred).^2) / denomPlot;
        end

        % ---------------------------------
        % Axis limits
        % ---------------------------------
        if isfinite(mn) && isfinite(mx) && mn < mx
            pad = 0.05 * (mx - mn);
            xlim([mn-pad mx+pad]);
            ylim([mn-pad mx+pad]);
        end

        xlabel('Observed Y');
        ylabel('Predicted Y');
        title(sprintf('Cross-validated predicted vs observed (N = %d)', numel(yObs)));
        grid on;
        set(gca,'FontSize',14,'LineWidth',1.2);
        axis square;

        % ---------------------------------
        % Annotation
        % ---------------------------------
        txt = {};
        if isfield(results,'Q2') && ~isempty(results.Q2)
            txt{end+1} = sprintf('Nested CV Q2 = %.3f', results.Q2);
        end
        if isfinite(rPlot)
            txt{end+1} = sprintf('Obs-Pred r = %.3f', rPlot);
        end
        if isfinite(R2plot)
            txt{end+1} = sprintf('Display R^2 = %.3f', R2plot);
        end
        if ~isempty(bOLS)
            txt{end+1} = sprintf('OLS slope = %.3f', bOLS(1));
        end
        if isfinite(robustSlope)
            txt{end+1} = sprintf('Robust slope = %.3f', robustSlope);
        end
        if isfield(results,'RMSE') && ~isempty(results.RMSE)
            txt{end+1} = sprintf('Mean RMSE = %.3f', results.RMSE);
        end
        if hasSubj
            txt{end+1} = sprintf('Unique subjects = %d', numel(unique(subjID)));
        end
        if hasRepeat
            txt{end+1} = sprintf('Repeats represented = %d', numel(unique(repID)));
        end

        if ~isempty(txt)
            text(0.05, 0.95, strjoin(txt, '\n'), 'Units', 'normalized', ...
                'VerticalAlignment', 'top', ...
                'BackgroundColor', 'w', ...
                'EdgeColor', [0.8 0.8 0.8], ...
                'Margin', 6);
        end

        legend('Location','bestoutside');
    end
end

% --------------------------
% 9. Permutation Q2 histogram (if available)
% --------------------------
if isfield(results,'allpermQ2') && ~isempty(results.allpermQ2)
    permQ2 = results.allpermQ2(:);
    permQ2 = permQ2(isfinite(permQ2));

    if ~isempty(permQ2)
        figure('Color','w','Units','pixels','Position',[100 100 900 700]);
        set(gcf,'Name','PLSR_permutation_Q2','NumberTitle','off');
        histogram(permQ2);
        hold on;
        if isfield(results,'Q2') && ~isempty(results.Q2)
            xline(results.Q2,'r-','LineWidth',2);
        end
        xlabel('Q2');
        ylabel('Frequency');
        title('Permutation Q2 distribution');
        grid on;
        set(gca,'FontSize',14,'LineWidth',1.2);

        txt = {};
        if isfield(results,'Q2') && ~isempty(results.Q2)
            txt{end+1} = sprintf('Observed Q2 = %.3f', results.Q2);
        end
        if isfield(results,'permQ2') && ~isempty(results.permQ2)
            txt{end+1} = sprintf('Mean perm Q2 = %.3f', results.permQ2);
        end
        if isfield(results,'permutation_p') && ~isempty(results.permutation_p)
            txt{end+1} = sprintf('p = %.3f', results.permutation_p);
        end
        if ~isempty(txt)
            text(0.05,0.95,strjoin(txt,'\n'),'Units','normalized', ...
                'VerticalAlignment','top','BackgroundColor','w');
        end
    end
end

% --------------------------
% 10. Bootstrap Q2 histogram (if available)
% --------------------------
if isfield(results,'allbootQ2') && ~isempty(results.allbootQ2)
    bootQ2 = results.allbootQ2(:);
    bootQ2 = bootQ2(isfinite(bootQ2));

    if ~isempty(bootQ2)
        figure('Color','w','Units','pixels','Position',[100 100 900 700]);
        set(gcf,'Name','PLSR_bootstrap_Q2','NumberTitle','off');
        histogram(bootQ2);
        hold on;
        if isfield(results,'Q2') && ~isempty(results.Q2)
            xline(results.Q2,'r-','LineWidth',2);
        end
        if isfield(results,'Q2_CI') && ~isempty(results.Q2_CI) && numel(results.Q2_CI)==2
            xline(results.Q2_CI(1),'k--','LineWidth',1.2);
            xline(results.Q2_CI(2),'k--','LineWidth',1.2);
        end
        xlabel('Q2');
        ylabel('Frequency');
        title('Bootstrap OOB Q2 distribution');
        grid on;
        set(gca,'FontSize',14,'LineWidth',1.2);

        txt = {};
        if isfield(results,'Q2') && ~isempty(results.Q2)
            txt{end+1} = sprintf('Observed Q2 = %.3f', results.Q2);
        end
        if isfield(results,'bootQ2') && ~isempty(results.bootQ2)
            txt{end+1} = sprintf('Mean boot Q2 = %.3f', results.bootQ2);
        end
        if isfield(results,'Q2_CI') && ~isempty(results.Q2_CI) && numel(results.Q2_CI)==2
            txt{end+1} = sprintf('95%% CI = [%.3f, %.3f]', results.Q2_CI(1), results.Q2_CI(2));
        end
        if ~isempty(txt)
            text(0.05,0.95,strjoin(txt,'\n'),'Units','normalized', ...
                'VerticalAlignment','top','BackgroundColor','w');
        end
    end
end

% --------------------------
% 11. Learning curve (if available)
% --------------------------
if isfield(results,'learningSizes') && isfield(results,'learningQ2') && ...
        ~isempty(results.learningSizes) && ~isempty(results.learningQ2)

    figure('Color','w','Units','pixels','Position',[100 100 1000 700]);
    set(gcf,'Name','PLSR_learning_curve','NumberTitle','off');
    plot(results.learningSizes(:), results.learningQ2(:), 'o-','LineWidth',1.5);
    xlabel('Sample size');
    ylabel('Q2');
    title('Learning curve');
    grid on;
    set(gca,'FontSize',14,'LineWidth',1.2);
end

end


% helper function

function img = orient_for_display(vol2d)
% Keep atlas, overlay, and underlay in the same display orientation
img = rot90(vol2d,-1);
end