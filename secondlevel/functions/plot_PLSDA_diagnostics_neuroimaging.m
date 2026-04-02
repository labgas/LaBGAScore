function ROI_table = plot_PLSDA_diagnostics_neuroimaging(results, XL, roiNames, atlasFile, varargin)
%
% Plot/diagnose PLS-DA neuroimaging results (ROI/parcellation).
%
% This function provides generic visualization and diagnostics for PLS-DA
% results produced by pipelines such as PLSDA_neuroimaging_pipeline /
% TSPO_PLSDA_pipeline. It creates:
%   - ROI table export (VIP + stabilityZ + robust flag)
%   - VIP vs stabilityZ scatter with robust ROI labeling
%   - NIfTI maps: VIP, stabilityZ, LV loadings (raw + thresholded)
%   - Multi-slice visualization (LV / VIP / stabilityZ) with ROI labels
%   - Optional: top-K mean weights table/bar plot (if results.meanFeatureWeight exists)
%   - Optional: feature stability stem plot (if results.featureStability exists)
%
% USAGE
%   ROI_table = plot_PLSDA_diagnostics_neuroimaging(results, XL, roiNames, atlasFile)
%   ROI_table = plot_PLSDA_diagnostics_neuroimaging(results, XL, roiNames, atlasFile, 'Name',Value,...)
%
% INPUTS
%   results struct
%        Must contain:
%          results.VIP        [p x 1] VIP scores from the final PLS model
%          results.stabilityZ [p x 1] stability statistic (meanBeta ./ sdBeta)
%        Optional fields used if present:
%          results.finalXLoadings    [p x finalLV] used if XL is empty
%          results.finalLV           scalar       clamps LV selection if present
%          results.meanFeatureWeight [p x 1] for top-K weight table/bar plot
%          results.selectionFrequency[p x 1] for top-K table (if available)
%          results.featureStability  [p x 1] for stability stem plot (if available)
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
%
% NAME–VALUE OPTIONS
%   'TopN'         (20)   Number of ROI-table rows printed to console (CSV exports all rows).
%   'TopK'         (20)   Number of top features for weight table/bar plot (if available).
%   'LV'           (2)    LV index to visualize (clamped to valid range).
%   'OutPrefix'    ('PLSDA')  Prefix for exported CSV/NIfTI files.
%   'VIP_thresh'   (1)    VIP threshold for robust contributor definition/labeling.
%   'stab_thresh'  (2)    stabilityZ threshold for robust contributor definition/labeling.
%   'MapPrctile'   (70)   Threshold LV map by abs percentile (keeps top 100-MapPrctile %).
%   'RelaxIfEmpty' (true) If no ROI passes thresholds, relax thresholds (75th percentiles)
%                         for visualization/labeling only.
%   'UnderlayFile' ('')   Optional structural underlay NIfTI for the multi-slice figure.
%                         Must be in the same voxel space/dimensions as atlasFile.
%
% OUTPUT
%   ROI_table table
%        Columns:
%          ROI, VIP, stabilityZ, RobustContributor
%        Also written to: <OutPrefix>_ROI_VIP_stability.csv
%
% FILES WRITTEN (requires SPM for NIfTI I/O)
%   <OutPrefix>_ROI_VIP_stability.csv
%   <OutPrefix>_VIP_map.nii
%   <OutPrefix>_stabilityZ_map.nii
%   <OutPrefix>_LV<LV>_map.nii
%   <OutPrefix>_LV<LV>_map_thresh.nii
%   (optional) <OutPrefix>_top<TopK>_weights.csv
%
% FIGURES CREATED (interactive)
%   1) VIP vs stabilityZ scatter (robust ROIs labeled)
%   2) Multi-slice 3×N panel: LV map / VIP map / stabilityZ map
%      with optional structural underlay and robust ROI labels
%   3) (optional) Top-K mean weight bar plot
%   4) (optional) Feature stability stem plot
%
% NOTES
%   - Robust contributors are defined as VIP>VIP_thresh AND |stabilityZ|>stab_thresh.
%   - Correct atlas/feature alignment is essential: atlas labels 1..p must match
%     feature order in X/results/XL.
%   - Figure sizes were enlarged relative to the original version to improve
%     readability when saved/exported.
%
% DEPENDENCIES
%   Requires SPM (spm_vol, spm_read_vols, spm_write_vol) for NIfTI mapping.
%
% See also: plsregress, perfcurve, spm_vol, spm_read_vols, spm_write_vol


% --------------------------
% Options
% --------------------------
p_in = inputParser;
addParameter(p_in,'TopN',20,@(x) isnumeric(x) && isscalar(x));
addParameter(p_in,'TopK',20,@(x) isnumeric(x) && isscalar(x));
addParameter(p_in,'LV',2,@(x) isnumeric(x) && isscalar(x));
addParameter(p_in,'OutPrefix','PLSDA',@(x) ischar(x) || isstring(x));
addParameter(p_in,'VIP_thresh',1,@(x) isnumeric(x) && isscalar(x));
addParameter(p_in,'stab_thresh',2,@(x) isnumeric(x) && isscalar(x));
addParameter(p_in,'MapPrctile',70,@(x) isnumeric(x) && isscalar(x));
addParameter(p_in,'RelaxIfEmpty',true,@(x) islogical(x) && isscalar(x));
addParameter(p_in,'UnderlayFile','',@(x) ischar(x) || isstring(x));
parse(p_in,varargin{:});

TopN         = p_in.Results.TopN;
TopK         = p_in.Results.TopK;
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
if ~isfield(results,'VIP') || isempty(results.VIP), error('results.VIP missing'); end
if ~isfield(results,'stabilityZ') || isempty(results.stabilityZ), error('results.stabilityZ missing'); end

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

% Robust LV selection/clamping
LV = max(1, min(LV, size(XL,2)));
% if isfield(results,'finalLV') && ~isempty(results.finalLV)
%     LV = max(1, min(results.finalLV, size(XL,2)));
% end

% --------------------------
% 1. ROI table
% --------------------------
isRobust = (results.VIP > VIP_thresh) & (abs(results.stabilityZ) > stab_thresh);

% Relax thresholds for visualization only if nothing passes
VIP_thresh_vis  = VIP_thresh;
stab_thresh_vis = stab_thresh;
if RelaxIfEmpty && ~any(isRobust)
    VIP_thresh_vis  = prctile(results.VIP,75);
    stab_thresh_vis = prctile(abs(results.stabilityZ),75);
    isRobust = (results.VIP > VIP_thresh_vis) & (abs(results.stabilityZ) > stab_thresh_vis);
end

ROI_table = table(roiNames, results.VIP(:), results.stabilityZ(:), isRobust(:), ...
    'VariableNames', {'ROI','VIP','stabilityZ','RobustContributor'});
ROI_table = sortrows(ROI_table,'stabilityZ','descend');

disp(ROI_table(1:min(TopN,height(ROI_table)),:))
writetable(ROI_table, sprintf('%s_ROI_VIP_stability.csv',OutPrefix));

% --------------------------
% 2. Scatter VIP vs stabilityZ
% --------------------------
figure('Color','w','Units','pixels','Position',[100 100 1000 800]); hold on;
scatter(results.VIP, results.stabilityZ, 60, 'filled');
xlabel('VIP','FontSize',14);
ylabel('stabilityZ','FontSize',14);
title('ROI importance vs stability (PLS-DA)','FontSize',16);
grid on;
set(gca,'FontSize',14,'LineWidth',1.2);
xline(VIP_thresh_vis,'r--','LineWidth',1.2);
yline(stab_thresh_vis,'r--','LineWidth',1.2);

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
V = spm_vol(atlasFile);
atlasData = spm_read_vols(V);

labels = unique(atlasData(:));
labels(labels==0 | isnan(labels)) = [];
if ~isempty(labels) && max(labels) < p
    warning('Atlas max label (%d) < number of features p (%d). Mapping may be wrong.', max(labels), p);
end

% Optional underlay
hasUnderlay = ~isempty(UnderlayFile) && exist(UnderlayFile,'file');
if hasUnderlay
    Vu = spm_vol(UnderlayFile);
    underlayData = spm_read_vols(Vu);

    if ~isequal(size(underlayData), size(atlasData))
        error('Underlay and atlas must have the same dimensions.');
    end
else
    underlayData = [];
end

% Prefer slices that contain atlas content
zHasData = squeeze(any(any(atlasData>0,1),2));
zList = find(zHasData);
if numel(zList) >= 5
    slice_idx = zList(round(linspace(1,numel(zList),5)));
else
    slice_idx = round(linspace(1,size(atlasData,3),5));
end

% --------------------------
% 4. VIP & stabilityZ maps
% --------------------------
VIPmap  = zeros(size(atlasData));
stabMap = zeros(size(atlasData));

for i = 1:p
    mask = atlasData == i;
    VIPmap(mask)  = results.VIP(i);
    stabMap(mask) = results.stabilityZ(i);
end

Vvip = V;  Vvip.fname  = sprintf('%s_VIP_map.nii',OutPrefix);
Vstab= V;  Vstab.fname = sprintf('%s_stabilityZ_map.nii',OutPrefix);
spm_write_vol(Vvip, VIPmap);
spm_write_vol(Vstab, stabMap);

% --------------------------
% 5. LV map (raw + thresholded)
% --------------------------
LVweights = XL(:,LV);
den = max(abs(LVweights));
if den==0 || isnan(den) || isinf(den), den = 1; end
LVweights = LVweights / den;

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
    if abs(val) < th, val = 0; end
    LVmap_thresh(mask) = val;
end
VlvT = V; VlvT.fname = sprintf('%s_LV%d_map_thresh.nii',OutPrefix,LV);
spm_write_vol(VlvT, LVmap_thresh);

% --------------------------
% 6. Multi-slice figure: LV / VIP / stabilityZ with optional underlay
%    (underlay kept as currently correct; overlays + labels flipped 180°)
% --------------------------

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
                text(axBase, median(xx), median(yy), roiNames{roiID}, ...
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

set(gcf,'Name','PLSDA_multislice','NumberTitle','off');

% --------------------------
% 7. Top-K weights + stability (if present)
% --------------------------
if isfield(results,'meanFeatureWeight') && ~isempty(results.meanFeatureWeight)
    mf = results.meanFeatureWeight(:);
    [~,idx] = sort(abs(mf),'descend');
    topK = min(TopK,p);

    if isfield(results,'selectionFrequency') && ~isempty(results.selectionFrequency)
        sf = results.selectionFrequency(:);
    else
        sf = nan(p,1);
    end

    weight_table = table(roiNames(idx(1:topK)), mf(idx(1:topK)), sf(idx(1:topK)), ...
        'VariableNames', {'feature label','mean weight','selection frequency'});
    disp(weight_table);
    writetable(weight_table, sprintf('%s_top%d_weights.csv',OutPrefix,topK));

    figure('Color','w','Units','pixels','Position',[100 100 1400 800]);
    bar(mf(idx(1:topK)));
    xticks(1:topK);
    xticklabels(roiNames(idx(1:topK)));
    set(gca,'TickLabelInterpreter','none');
    xtickangle(45);
    ylabel('Mean Weight','FontSize',14);
    title(sprintf('Top %d Features (mean weight)',topK),'FontSize',16);
    set(gca,'FontSize',14,'LineWidth',1.2);
    set(gcf,'Name','Top20_weights','NumberTitle','off');
end

if isfield(results,'featureStability') && ~isempty(results.featureStability)
    figure('Color','w','Units','pixels','Position',[100 100 1400 700]);
    stem(results.featureStability(:),'LineWidth',1.2);
    xlabel('Feature Index','FontSize',14);
    ylabel('Stability Proportion','FontSize',14);
    title('Feature Selection Stability','FontSize',16);
    set(gca,'FontSize',14,'LineWidth',1.2);
    set(gcf,'Name','Feature_stability','NumberTitle','off');
end

end


% helper function

function img = orient_for_display(vol2d)
% Keep atlas, overlay, and underlay in the same display orientation
img = rot90(vol2d,-1);
end