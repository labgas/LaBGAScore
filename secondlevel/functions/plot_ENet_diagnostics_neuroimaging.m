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
addParameter(p_input,'FreqThresh',0.5,@(x) isnumeric(x) && isscalar(x));
addParameter(p_input,'WeightThresh',0,@(x) isnumeric(x) && isscalar(x));
addParameter(p_input,'MapPrctile',70,@(x) isnumeric(x) && isscalar(x));
addParameter(p_input,'DoPostSelection',true,@(x) islogical(x) && isscalar(x));
addParameter(p_input,'OutPrefix','ENet',@(x) ischar(x) || isstring(x));
addParameter(p_input,'RelaxIfEmpty',true,@(x) islogical(x) && isscalar(x));
addParameter(p_input,'UnderlayFile','',@(x) ischar(x) || isstring(x));
parse(p_input,varargin{:});

TopN          = p_input.Results.TopN;
FreqThresh    = p_input.Results.FreqThresh;
WeightThresh  = p_input.Results.WeightThresh;
MapPrctile    = p_input.Results.MapPrctile;
DoPostSel     = p_input.Results.DoPostSelection;
OutPrefix     = char(p_input.Results.OutPrefix);
RelaxIfEmpty  = p_input.Results.RelaxIfEmpty;
UnderlayFile  = char(p_input.Results.UnderlayFile);

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
roiNames = cellstr(roiNames(:));

if ~isfield(results,'meanFeatureWeight')
    error('results.meanFeatureWeight missing. Ensure ENet pipeline stores meanFeatureWeight.');
end
if ~isfield(results,'featureWeights')
    error('results.featureWeights missing. Ensure ENet pipeline stores featureWeights.');
end
if ~isfield(results,'selectionFrequency')
    warning('results.selectionFrequency missing: computing TopK frequency as zeros.');
    results.selectionFrequency = zeros(p,1);
end
if ~isfield(results,'featureStability')
    results.featureStability = mean(abs(results.featureWeights)>0,2);
end

if isfield(results,'selectionTopK') && ~isempty(results.selectionTopK)
    TopK = results.selectionTopK;
else
    % fallback for backward compatibility
    p = length(results.meanFeatureWeight);
    TopK = min(20, max(3, ceil(0.25 * p)));
    TopK = min(TopK, p-1);
    TopK = max(1, TopK);
end

% Convert Y to 0/1 numeric
if iscell(Y) || isstring(Y) || iscategorical(Y)
    Y = grp2idx(Y);
end
yNum = double(Y(:)==max(Y));

meanW = results.meanFeatureWeight(:);
absW  = abs(meanW);
selFreq  = results.selectionFrequency(:);
stabProp = results.featureStability(:);

% Robustness criterion
isRobust = (absW > WeightThresh) & (selFreq >= FreqThresh);

% Relax thresholds for visualization/labels only if none robust
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
figure('Color','w','Units','pixels','Position',[100 100 1000 800]); hold on;
set(gcf,'Name','ENet_meanWeight_vs_selectionFrequency','NumberTitle','off');

scatter(absW, selFreq, 60, 'b', 'filled');
xlabel('|Mean ENet weight|','FontSize',14);
ylabel('Selection frequency (TopK fraction)','FontSize',14);
title('ROI importance and robustness (Elastic Net)','FontSize',16);
grid on;
set(gca,'FontSize',14,'LineWidth',1.2);

try
    hy = yline(FreqThresh_vis,'r--','LineWidth',1.2);
    hx = xline(WeightThresh_vis,'r--','LineWidth',1.2);
    hy.Label = sprintf('freq >= %.2f', FreqThresh_vis);
    hx.Label = sprintf('|w| > %.2f', WeightThresh_vis);
catch
    yline(FreqThresh_vis,'r--',sprintf('freq >= %.2f', FreqThresh_vis),'LineWidth',1.2);
    xline(WeightThresh_vis,'r--',sprintf('|w| > %.2f', WeightThresh_vis),'LineWidth',1.2);
end

robustIdx = find(isRobust);
for ii = 1:numel(robustIdx)
    idx = robustIdx(ii);
    text(absW(idx)+0.02*max(absW), selFreq(idx), roiNames{idx}, ...
        'FontSize',11, 'Interpreter','none');
end

% --------------------------
% 3. Convert meanWeight and selectionFrequency to NIfTI
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
% 4. Multi-slice figure: meanWeight / selectionFrequency / featureStability
% --------------------------

% Recompute slice selection and avoid extreme inferior/superior slices
zHasData = squeeze(any(any(atlasData>0,1),2));
zList = find(zHasData);

if numel(zList) >= 5
    lo = max(1, round(0.20 * numel(zList)));
    hi = min(numel(zList), round(0.80 * numel(zList)));
    keepIdx = zList(lo:hi);

    if numel(keepIdx) >= 5
        slice_idx = keepIdx(round(linspace(1, numel(keepIdx), 5)));
    else
        slice_idx = zList(round(linspace(1, numel(zList), 5)));
    end
else
    slice_idx = round(linspace(1,size(wMap,3),5));
end

% Display maps
wDisp = wMap;
fDisp = fMap; 
fDisp(fDisp < FreqThresh_vis) = 0;
sDisp = sMap;

figure('Color','k','Units','pixels','Position',[50 50 2200 1200]);
set(gcf,'Name','ENet_multislice','NumberTitle','off');

nSlices = length(slice_idx);
baseAxes = gobjects(3,nSlices);

for s = 1:nSlices
    z = slice_idx(s);

    % ---------- Row 1: meanWeight ----------
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
    ov = rot90(orient_for_display(wDisp(:,:,z)), 2);
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

    % ---------- Row 2: selectionFrequency ----------
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
    ov = rot90(orient_for_display(fDisp(:,:,z)), 2);
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

    % ---------- Row 3: featureStability ----------
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
    ov = rot90(orient_for_display(sDisp(:,:,z)), 2);
    hOv = imagesc(ov,'Parent',axOverlay);
    axis(axOverlay,'image');
    axis(axOverlay,'off');
    colormap(axOverlay,jet);

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
end

% ROI labels
for s = 1:nSlices
    z = slice_idx(s);

    for row = 1:3
        axBase = baseAxes(row,s);
        hold(axBase,'on');

        for ii = 1:length(robustIdx)
            roiID = robustIdx(ii);

            mask = rot90(orient_for_display(atlasData(:,:,z) == roiID), 2);
            [yy,xx] = find(mask);

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

% Column headers: z in mm
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
rowNames = {'MeanW','Freq','Stab'};
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

% --------------------------
% 5. Top-K weight table & bar plot
% --------------------------
[~,idx] = sort(absW,'descend');

weight_table = table(roiNames(idx(1:TopK)), ...
    meanW(idx(1:TopK)), ...
    selFreq(idx(1:TopK)), ...
    stabProp(idx(1:TopK)), ...
    'VariableNames', {'feature label','mean weight','selection frequency','feature stability'});

disp(weight_table);

figure('Color','w','Units','pixels','Position',[100 100 1400 800]);
set(gcf,'Name',sprintf('ENet_top%d_weights',TopK),'NumberTitle','off');

bar(meanW(idx(1:TopK)));
xticks(1:TopK);
xticklabels(roiNames(idx(1:TopK)));
set(gca,'TickLabelInterpreter','none');
xtickangle(45);
ylabel('Mean Weight','FontSize',14);
title(sprintf('Top %d Stable Brain Features (Elastic Net)',TopK),'FontSize',16);
set(gca,'FontSize',14,'LineWidth',1.2);

% --------------------------
% 6. Feature selection stability plot
% --------------------------
figure('Color','w','Units','pixels','Position',[100 100 1400 700]);
set(gcf,'Name','ENet_feature_stability','NumberTitle','off');

stem(stabProp,'LineWidth',1.2);
xlabel('Feature Index','FontSize',14);
ylabel('Stability proportion (non-zero across CV runs)','FontSize',14);
title('Feature Selection Stability','FontSize',16);
set(gca,'FontSize',14,'LineWidth',1.2);

% --------------------------
% 7. Alpha/Lambda selection histograms (if available)
% --------------------------
if isfield(results,'selectedAlpha') && ~isempty(results.selectedAlpha)
    figure('Color','w','Units','pixels','Position',[100 100 900 700]);
    set(gcf,'Name','ENet_alpha_selection_histogram','NumberTitle','off');

    histogram(results.selectedAlpha(:));
    xlabel('Selected alpha','FontSize',14);
    ylabel('Frequency','FontSize',14);
    title('Alpha selection across CV folds','FontSize',16);
    set(gca,'FontSize',14,'LineWidth',1.2);
end

if isfield(results,'selectedLambda') && ~isempty(results.selectedLambda)
    figure('Color','w','Units','pixels','Position',[100 100 900 700]);
    set(gcf,'Name','ENet_lambda_selection_histogram','NumberTitle','off');

    histogram(log10(results.selectedLambda(:)));
    xlabel('log10(selected lambda)','FontSize',14);
    ylabel('Frequency','FontSize',14);
    title('Lambda selection across CV folds','FontSize',16);
    set(gca,'FontSize',14,'LineWidth',1.2);
end

% --------------------------
% 8. Post-selection inference: betas & p-values for top non-zero regions
% --------------------------
if DoPostSel
    nz = find(meanW ~= 0);
    if ~isempty(nz)
        [~,ord] = sort(absW(nz),'descend');
        keep = nz(ord(1:min(TopK,numel(nz))));

        Xsel = X(:,keep);

        try
            % Newer MATLAB versions
            mdl = fitglm(Xsel, yNum, ...
                'Distribution','binomial', ...
                'Link','logit', ...
                'LikelihoodPenalty','jeffreys-prior');
        catch ME
            if contains(ME.message,'LikelihoodPenalty') || contains(ME.identifier,'parseArgs')
                warning('LikelihoodPenalty not supported in this MATLAB version; falling back to standard logistic regression.');
                mdl = fitglm(Xsel, yNum, ...
                    'Distribution','binomial', ...
                    'Link','logit');
            else
                rethrow(ME);
            end
        end

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

function img = orient_for_display(vol2d)
img = rot90(vol2d, -1);
end