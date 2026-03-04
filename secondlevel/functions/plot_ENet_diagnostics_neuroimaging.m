function ROI_table = plot_ENet_diagnostics_neuroimaging(results, X, Y, roiNames, atlasFile, varargin)

% Plotting + diagnostics for Elastic Net neuroimaging pipelines.
%
% This function visualizes and exports diagnostics from Elastic Net pipelines
% (e.g., TSPO_ENet_pipeline / ENet_neuroimaging_pipeline). It mirrors the
% architecture of the TSPO-specific ENet plotting script while remaining generic
% to any neuroimaging feature matrix where columns correspond to atlas ROI labels
% (1..p). Outputs include:
%   - ROI table export (mean weights + selection frequency + stability + robust flag)
%   - |mean weight| vs selection frequency scatter with robust ROI labels
%   - NIfTI maps: mean weights, selection frequency, feature stability (raw + thresholded)
%   - Multi-slice axial visualization (weights / selection freq / stability)
%   - Top-K weights bar plot + feature stability stem plot
%   - Alpha/Lambda selection histograms (if present)
%   - Optional post-selection logistic refit (betas/SE/OR/p-values) for top non-zero regions
%
% USAGE
%   ROI_table = plot_ENet_diagnostics_neuroimaging(results, X, Y, roiNames, atlasFile)
%   ROI_table = plot_ENet_diagnostics_neuroimaging(results, X, Y, roiNames, atlasFile, 'Name',Value,...)
%
% INPUTS
%   results   struct
%            Output struct from an ENet pipeline. Must contain:
%              results.meanFeatureWeight [p x 1] mean ENet beta across CV runs
%              results.featureWeights    [p x R] beta matrix across runs (R = repeats*folds)
%            Typically also contains:
%              results.selectionFrequency [p x 1] top-K selection frequency across runs
%              results.featureStability   [p x 1] proportion of runs where beta ~= 0
%              results.selectedAlpha      [nRepeats x outerK] (optional; histogram)
%              results.selectedLambda     [nRepeats x outerK] (optional; histogram)
%
%   X        [n x p] numeric
%            Feature matrix used in the pipeline. Required to (a) validate p,
%            and (b) run optional post-selection logistic refit.
%
%   Y        [n x 1] labels
%            Binary labels (numeric/logical/categorical/string/cellstr).
%            Converted internally using:
%              yNum = double(Y == max(Y))
%
%   roiNames [p x 1] cellstr/string (optional)
%            ROI names aligned to feature order. If empty/missing, names are
%            auto-generated as ROI_001..ROI_p.
%
%   atlasFile char/string
%            Path to an atlas NIfTI where ROI labels are integers 1..p (0 background).
%            Used to write NIfTI maps for ENet-derived quantities.
%
% NAME-VALUE OPTIONS
%   'TopN'            (default 20)  Number of top rows printed/exported in ROI table.
%   'TopK'            (default 20)  Number of top features in bar plot/table and post-selection refit.
%   'FreqThresh'      (default 0.5) Robustness threshold on selectionFrequency.
%   'WeightThresh'    (default 0)   Robustness threshold on |meanFeatureWeight|.
%   'MapPrctile'      (default 70)  Threshold mean-weight map by abs percentile.
%   'DoPostSelection' (default true) Run post-selection logistic refit on top non-zero features.
%   'OutPrefix'       (default 'ENet') Prefix for exported files.
%   'RelaxIfEmpty'    (default true) If no ROI passes thresholds, relax thresholds
%                                  for visualization/labels ONLY (75th percentiles).
%
% OUTPUT
%   ROI_table table
%            Table with columns:
%              ROI, meanWeight, absMeanWeight, selectionFrequency, featureStability, RobustContributor
%            Sorted descending by absMeanWeight and exported to CSV.
%
% FILES WRITTEN (SPM REQUIRED)
%   {OutPrefix}_ROI_weights_stability.csv
%   {OutPrefix}_meanWeight_map.nii
%   {OutPrefix}_selectionFreq_map.nii
%   {OutPrefix}_featureStability_map.nii
%   {OutPrefix}_meanWeight_map_thresh.nii
%   {OutPrefix}_postselection_refit_table.csv (if DoPostSelection and any non-zero weights)
%
% FIGURES GENERATED
%   1) |meanWeight| vs selectionFrequency scatter (robust ROIs labeled)
%   2) Multi-slice axial figure (rows: meanWeight / selectionFreq / featureStability)
%   3) Top-K meanWeight bar plot
%   4) Feature stability stem plot
%   5) Optional: selectedAlpha histogram (if results.selectedAlpha present)
%   6) Optional: log10(selectedLambda) histogram (if results.selectedLambda present)
%
% INTERPRETATION NOTES
%   - meanFeatureWeight summarizes the direction/magnitude of effects across CV runs.
%   - selectionFrequency is a robustness metric: how often the feature is among top-K |beta|.
%   - featureStability captures sparsity: proportion of runs where beta ~= 0.
%   - Post-selection refit is exploratory: it conditions on a data-driven selection step.
%     Jeffreys prior penalty is used to mitigate separation, but p-values should be
%     interpreted cautiously (not confirmatory).
%
% IMPLEMENTATION NOTES
%   - Atlas sanity check warns if max atlas label < p (mapping may be wrong).
%   - Slice selection prefers slices that contain atlas content.
%   - RelaxIfEmpty relaxes thresholds only for visualization/labels to avoid empty maps.
%
% DEPENDENCIES
%   Requires SPM12 on the MATLAB path for NIfTI I/O:
%     spm_vol, spm_read_vols, spm_write_vol
%   Requires Statistics and Machine Learning Toolbox:
%     table, fitglm
%
% See also: lasso, fitglm, perfcurve, spm_write_vol

% --------------------------
% Parse options
% --------------------------
p_input = inputParser;
addParameter(p_input,'TopN',20,@(x) isnumeric(x) && isscalar(x));
addParameter(p_input,'TopK',20,@(x) isnumeric(x) && isscalar(x));
addParameter(p_input,'FreqThresh',0.5,@(x) isnumeric(x) && isscalar(x));     % robustness threshold for selectionFrequency
addParameter(p_input,'WeightThresh',0,@(x) isnumeric(x) && isscalar(x));     % robustness threshold for |meanWeight|
addParameter(p_input,'MapPrctile',70,@(x) isnumeric(x) && isscalar(x));      % threshold maps at top (100-MapPrctile)% abs weights
addParameter(p_input,'DoPostSelection',true,@(x) islogical(x) && isscalar(x));
addParameter(p_input,'OutPrefix','ENet',@(x) ischar(x) || isstring(x));
addParameter(p_input,'RelaxIfEmpty',true,@(x) islogical(x) && isscalar(x)); % NEW: relax thresholds if none robust
parse(p_input,varargin{:});

TopN          = p_input.Results.TopN;
TopK          = p_input.Results.TopK;
FreqThresh    = p_input.Results.FreqThresh;
WeightThresh  = p_input.Results.WeightThresh;
MapPrctile    = p_input.Results.MapPrctile;
DoPostSel     = p_input.Results.DoPostSelection;
OutPrefix     = char(p_input.Results.OutPrefix);
RelaxIfEmpty  = p_input.Results.RelaxIfEmpty;

% --------------------------
% Basic checks
% --------------------------
if nargin < 2 || isempty(X)
    error('ENet_plot_neuroimaging requires X (n x p) as used in ENet pipeline.');
end
p = size(X,2);

if nargin < 3 || isempty(Y)
    error('ENet_plot_neuroimaging requires Y (n x 1) outcome labels as used in ENet pipeline.');
end

if ~exist('roiNames','var') || isempty(roiNames)
    roiNames = arrayfun(@(i) sprintf('ROI_%03d',i), 1:p, 'UniformOutput', false);
end

if ~isfield(results,'meanFeatureWeight')
    error('results.meanFeatureWeight missing. Ensure ENet pipeline stores meanFeatureWeight.');
end
if ~isfield(results,'featureWeights')
    error('results.featureWeights missing. Ensure ENet pipeline stores featureWeights.');
end
if ~isfield(results,'selectionFrequency')
    % fallback if not present
    warning('results.selectionFrequency missing: computing TopK frequency as zeros.');
    results.selectionFrequency = zeros(p,1);
end
if ~isfield(results,'featureStability')
    % fallback if not present
    results.featureStability = mean(abs(results.featureWeights)>0,2);
end

% Convert Y to 0/1 numeric (same logic as your pipelines)
if iscell(Y) || isstring(Y) || iscategorical(Y)
    Y = grp2idx(Y);
end
yNum = double(Y(:)==max(Y));

meanW = results.meanFeatureWeight(:);
absW  = abs(meanW);

% "Stability" analogue: use selectionFrequency and featureStability (nonzero proportion)
selFreq  = results.selectionFrequency(:);
stabProp = results.featureStability(:);

% Robustness criterion
isRobust = (absW > WeightThresh) & (selFreq >= FreqThresh);

% NEW: relax thresholds for visualization/labels only if nothing is robust
FreqThresh_vis   = FreqThresh;
WeightThresh_vis = WeightThresh;

if RelaxIfEmpty && ~any(isRobust)
    WeightThresh_vis = prctile(absW,75);
    FreqThresh_vis   = prctile(selFreq,75);
    isRobust = (absW > WeightThresh_vis) & (selFreq >= FreqThresh_vis);
end

% --------------------------
% 1. ROI reliability table
% --------------------------
ROI_table = table(roiNames(:), meanW, absW, selFreq, stabProp, isRobust, ...
    'VariableNames', {'ROI','meanWeight','absMeanWeight','selectionFrequency','featureStability','RobustContributor'});

ROI_table = sortrows(ROI_table, 'absMeanWeight', 'descend');

disp(ROI_table(1:min(TopN,height(ROI_table)),:))

csvName = sprintf('%s_ROI_weights_stability.csv',OutPrefix);
writetable(ROI_table, csvName);
disp(['Table exported: ' csvName]);

% --------------------------
% 2. Scatter: |meanWeight| vs selectionFrequency
% --------------------------
figure('Color','w'); hold on;
scatter(absW, selFreq, 50, 'b', 'filled');
xlabel('|Mean ENet weight|');
ylabel('Selection frequency (TopK fraction)');
title('ROI importance and robustness (Elastic Net)');
grid on;

yline(FreqThresh_vis,'r--',sprintf('freq \\ge %.2f',FreqThresh_vis));
xline(WeightThresh_vis,'r--',sprintf('|w| > %.2f',WeightThresh_vis));

robustIdx = find(isRobust);
for ii = 1:numel(robustIdx)
    idx = robustIdx(ii);
    text(absW(idx)+0.02*max(absW), selFreq(idx), roiNames{idx}, 'FontSize',10);
end

% --------------------------
% 3. Convert meanWeight and selectionFrequency to NIfTI
% --------------------------
V = spm_vol(atlasFile);
atlasData = spm_read_vols(V);

% NEW: atlas label sanity warning
labels = unique(atlasData(:));
labels(labels==0 | isnan(labels)) = [];
if ~isempty(labels) && max(labels) < p
    warning('Atlas max label (%d) < number of features p (%d). Mapping may be wrong.', max(labels), p);
end

wMap  = zeros(size(atlasData));
fMap  = zeros(size(atlasData));
sMap  = zeros(size(atlasData));

for i = 1:p
    mask = atlasData==i;
    wMap(mask) = meanW(i);
    fMap(mask) = selFreq(i);
    sMap(mask) = stabProp(i);
end

Vw = V; Vw.fname = sprintf('%s_meanWeight_map.nii',OutPrefix);
spm_write_vol(Vw,wMap); disp(['Mean weight NIfTI saved: ' Vw.fname]);

Vf = V; Vf.fname = sprintf('%s_selectionFreq_map.nii',OutPrefix);
spm_write_vol(Vf,fMap); disp(['Selection frequency NIfTI saved: ' Vf.fname]);

Vs = V; Vs.fname = sprintf('%s_featureStability_map.nii',OutPrefix);
spm_write_vol(Vs,sMap); disp(['Feature stability NIfTI saved: ' Vs.fname]);

% Thresholded weight map (~top (100-MapPrctile)%)
threshW = prctile(abs(meanW),MapPrctile);
wMap_thresh = zeros(size(atlasData));
for i = 1:p
    mask = atlasData==i;
    val = meanW(i);
    if abs(val) < threshW, val = 0; end
    wMap_thresh(mask) = val;
end
Vwt = V; Vwt.fname = sprintf('%s_meanWeight_map_thresh.nii',OutPrefix);
spm_write_vol(Vwt,wMap_thresh);
disp(['Thresholded mean weight NIfTI saved: ' Vwt.fname]);

% --------------------------
% 4. Multi-panel figure: weight / selectionFrequency / featureStability
% --------------------------

% NEW: choose slices that actually contain atlas content
zHasData = squeeze(any(any(atlasData>0,1),2));
zList = find(zHasData);
if numel(zList) >= 5
    slice_idx = zList(round(linspace(1,numel(zList),5)));
else
    slice_idx = round(linspace(1,size(wMap,3),5));
end

% Apply thresholds for display
wDisp = wMap; % optionally threshold later
fDisp = fMap; fDisp(fDisp < FreqThresh_vis) = 0;
sDisp = sMap;

figure('Color','w','Position',[50 50 1200 450]);
for s = 1:length(slice_idx)
    z = slice_idx(s);

    subplot(3,length(slice_idx),s);
    imagesc(wDisp(:,:,z)'); axis image; axis off;
    title(sprintf('MeanW Z=%d',z));

    subplot(3,length(slice_idx),s+length(slice_idx));
    imagesc(fDisp(:,:,z)'); axis image; axis off;
    title(sprintf('Freq Z=%d',z));

    subplot(3,length(slice_idx),s+2*length(slice_idx));
    imagesc(sDisp(:,:,z)'); axis image; axis off;
    title(sprintf('Stab Z=%d',z));

    % ROI labels for robust contributors
    for ii = 1:length(robustIdx)
        roiID = robustIdx(ii);
        mask = atlasData(:,:,z)==roiID;
        [y,x] = find(mask);
        if ~isempty(x)
            text(median(x), median(y), roiNames{roiID}, 'Color','w', ...
                'FontSize',8,'FontWeight','bold','HorizontalAlignment','center');
        end
    end
end
colormap('jet');
try
    sgtitle('Elastic Net: mean weights, selection frequency, feature stability (robust ROIs labeled)');
catch
    subtitle('Elastic Net: mean weights, selection frequency, feature stability (robust ROIs labeled)');
end
colorbar;

% --------------------------
% 5. Top-K weight table & bar plot
% --------------------------
[~,idx] = sort(absW,'descend');
topK = min(TopK,p);

weight_table = table(roiNames(idx(1:topK)), ...
    meanW(idx(1:topK)), ...
    selFreq(idx(1:topK)), ...
    stabProp(idx(1:topK)), ...
    'VariableNames', {'feature label','mean weight','selection frequency','feature stability'});

disp(weight_table);

figure('Color','w');
bar(meanW(idx(1:topK)));
xticks(1:topK);
xticklabels(roiNames(idx(1:topK)));
xtickangle(45);
ylabel('Mean Weight');
title(sprintf('Top %d Stable Brain Features (Elastic Net)',topK));

% --------------------------
% 6. Feature selection stability plot
% --------------------------
figure('Color','w');
stem(stabProp);
xlabel('Feature Index');
ylabel('Stability proportion (non-zero across CV runs)');
title('Feature Selection Stability');

% --------------------------
% 7. Alpha/Lambda selection histograms (if available)
% --------------------------
if isfield(results,'selectedAlpha')
    figure('Color','w');
    histogram(results.selectedAlpha(:));
    xlabel('Selected alpha'); ylabel('Frequency');
    title('Alpha selection across CV folds');
end
if isfield(results,'selectedLambda')
    figure('Color','w');
    histogram(log10(results.selectedLambda(:)));
    xlabel('log10(selected lambda)'); ylabel('Frequency');
    title('Lambda selection across CV folds');
end

% --------------------------
% 8. Post-selection inference: betas & p-values for top non-zero regions
% --------------------------
if DoPostSel
    nz = find(meanW ~= 0);
    if ~isempty(nz)
        [~,ord] = sort(absW(nz),'descend');
        keep = nz(ord(1:min(topK,numel(nz))));

        Xsel = X(:,keep);

        % Version-robust logistic regression + Firth/Jeffreys penalty to handle separation
        mdl = fitglm(Xsel, yNum, ...
            'Distribution','binomial', ...
            'Link','logit', ...
            'LikelihoodPenalty','jeffreys-prior');

        c = mdl.Coefficients;
        beta = c.Estimate(2:end);
        se   = c.SE(2:end);
        pval = c.pValue(2:end);
        OR   = exp(beta);

        postSel_table = table(roiNames(keep), beta, se, OR, pval, ...
            meanW(keep), selFreq(keep), stabProp(keep), ...
            'VariableNames', {'ROI','Beta_refit','SE','OddsRatio','pValue', ...
                              'ENet_meanWeight','selectionFrequency','featureStability'});

        postSel_table = sortrows(postSel_table,'pValue','ascend');

        disp('Post-selection refit (logistic regression) on top non-zero ENet regions:');
        disp(postSel_table);

        writetable(postSel_table, sprintf('%s_postselection_refit_table.csv',OutPrefix));
        disp(['Post-selection table exported: ' sprintf('%s_postselection_refit_table.csv',OutPrefix)]);
    else
        disp('No non-zero mean weights found; skipping post-selection inference table.');
    end
end

end