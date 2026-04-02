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
%          results.finalXLoadings   [p x finalLV] used if XL is empty
%          results.finalLV          scalar       clamps LV selection if present
%          results.meanFeatureWeight [p x 1] for top-K weight table/bar plot
%          results.selectionFrequency [p x 1] for top-K table (if available)
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
%   'TopN'        (20)   Number of ROI-table rows printed to console (CSV exports all rows).
%   'TopK'        (20)   Number of top features for weight table/bar plot (if available).
%   'LV'          (2)    LV index to visualize (clamped to valid range).
%   'OutPrefix'   ('PLSDA')  Prefix for exported CSV/NIfTI files.
%   'VIP_thresh'  (1)    VIP threshold for robust contributor definition/labeling.
%   'stab_thresh' (2)    stabilityZ threshold for robust contributor definition/labeling.
%   'MapPrctile'  (70)   Threshold LV map by abs percentile (keeps top 100-MapPrctile %).
%   'RelaxIfEmpty'(true) If no ROI passes thresholds, relax thresholds (75th percentiles)
%                        for visualization/labeling only.
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
%   2) Multi-slice 3×N panel: LV map / VIP map / stabilityZ map (robust ROIs labeled)
%   3) (optional) Top-K mean weight bar plot
%   4) (optional) Feature stability stem plot
%
% NOTES
%   - Robust contributors are defined as VIP>VIP_thresh AND |stabilityZ|>stab_thresh.
%   - Correct atlas/feature alignment is essential: atlas labels 1..p must match
%     feature order in X/results/XL.
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
parse(p_in,varargin{:});

TopN        = p_in.Results.TopN;
TopK        = p_in.Results.TopK;
LV          = p_in.Results.LV;
OutPrefix   = char(p_in.Results.OutPrefix);
VIP_thresh  = p_in.Results.VIP_thresh;
stab_thresh = p_in.Results.stab_thresh;
MapPrctile  = p_in.Results.MapPrctile;
RelaxIfEmpty= p_in.Results.RelaxIfEmpty;

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
roiNames = roiNames(:);

% Robust LV selection/clamping
LV = max(1, min(LV, size(XL,2)));
% if isfield(results,'finalLV') && ~isempty(results.finalLV)
%     LV = max(1, min(results.finalLV, size(XL,2)));
% end

% --------------------------
% 1. ROI table
% --------------------------
isRobust = (results.VIP>VIP_thresh) & (abs(results.stabilityZ)>stab_thresh);

% Relax thresholds for visualization only if nothing passes
VIP_thresh_vis = VIP_thresh;
stab_thresh_vis = stab_thresh;
if RelaxIfEmpty && ~any(isRobust)
    VIP_thresh_vis  = prctile(results.VIP,75);
    stab_thresh_vis = prctile(abs(results.stabilityZ),75);
    isRobust = (results.VIP>VIP_thresh_vis) & (abs(results.stabilityZ)>stab_thresh_vis);
end

ROI_table = table(roiNames, results.VIP(:), results.stabilityZ(:), isRobust(:), ...
    'VariableNames', {'ROI','VIP','stabilityZ','RobustContributor'});
ROI_table = sortrows(ROI_table,'stabilityZ','descend');

disp(ROI_table(1:min(TopN,height(ROI_table)),:))
writetable(ROI_table, sprintf('%s_ROI_VIP_stability.csv',OutPrefix));

% --------------------------
% 2. Scatter VIP vs stabilityZ
% --------------------------
figure('Color','w'); hold on;
scatter(results.VIP, results.stabilityZ, 50, 'filled');
xlabel('VIP'); ylabel('stabilityZ');
title('ROI importance vs stability (PLS-DA)'); grid on;
xline(VIP_thresh_vis,'r--'); yline(stab_thresh_vis,'r--');

robustIdx = find(isRobust);
for ii=1:length(robustIdx)
    idx = robustIdx(ii);
    text(results.VIP(idx)+0.02, results.stabilityZ(idx), roiNames{idx}, 'FontSize',9);
end

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

for i=1:p
    mask = atlasData==i;
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
for i=1:p
    mask = atlasData==i;
    LVmap(mask) = LVweights(i);
end

Vlv = V; Vlv.fname = sprintf('%s_LV%d_map.nii',OutPrefix,LV);
spm_write_vol(Vlv, LVmap);

th = prctile(abs(LVweights), MapPrctile);
LVmap_thresh = zeros(size(atlasData));
for i=1:p
    mask = atlasData==i;
    val = LVweights(i);
    if abs(val) < th, val = 0; end
    LVmap_thresh(mask) = val;
end
VlvT = V; VlvT.fname = sprintf('%s_LV%d_map_thresh.nii',OutPrefix,LV);
spm_write_vol(VlvT, LVmap_thresh);

% --------------------------
% 6. Multi-slice figure: LV / VIP / stabilityZ
% --------------------------
VIPmap_vis = VIPmap;  VIPmap_vis(VIPmap_vis<VIP_thresh_vis)=0;
stabMap_vis= stabMap; stabMap_vis(abs(stabMap_vis)<stab_thresh_vis)=0;

figure('Color','w','Position',[50 50 1200 450]);
for s=1:length(slice_idx)
    z = slice_idx(s);

    subplot(3,length(slice_idx),s)
    imagesc(LVmap(:,:,z)'); axis image off
    title(sprintf('LV%d Z=%d',LV,z))

    subplot(3,length(slice_idx),s+length(slice_idx))
    imagesc(VIPmap_vis(:,:,z)'); axis image off
    title(sprintf('VIP Z=%d',z))

    subplot(3,length(slice_idx),s+2*length(slice_idx))
    imagesc(stabMap_vis(:,:,z)'); axis image off
    title(sprintf('stabZ Z=%d',z))

    for ii=1:length(robustIdx)
        roiID = robustIdx(ii);
        mask = atlasData(:,:,z)==roiID;
        [y,x] = find(mask);
        if ~isempty(x)
            text(median(x),median(y),roiNames{roiID},'Color','w', ...
                'FontSize',8,'FontWeight','bold','HorizontalAlignment','center');
        end
    end
end
colormap('jet'); colorbar;
try
    sgtitle(sprintf('PLS-DA LV%d / VIP / stabilityZ (robust ROIs labeled)',LV));
catch
    subtitle(sprintf('PLS-DA LV%d / VIP / stabilityZ (robust ROIs labeled)',LV));
end

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

    figure('Color','w');
    bar(mf(idx(1:topK)));
    xticks(1:topK); xticklabels(roiNames(idx(1:topK)));
    xtickangle(45);
    ylabel('Mean Weight');
    title(sprintf('Top %d Features (mean weight)',topK));
end

if isfield(results,'featureStability') && ~isempty(results.featureStability)
    figure('Color','w');
    stem(results.featureStability(:));
    xlabel('Feature Index');
    ylabel('Stability Proportion');
    title('Feature Selection Stability');
end

end