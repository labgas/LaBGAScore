%% LaBGAScore_pdm_regenerate_reports.m
%
% *USAGE*
%
% Regenerates c2a-style reporting (overview montage, labelled region table,
% regioncenters montage) for PDM analyses that have already been run, reading
% the written PDM*.nii rather than refitting. Each bootstrap takes over an
% hour, so refitting purely to improve the figures is not worth it.
%
% *NOTES*
%
% The number of suprathreshold voxels is counted from the image, not read back
% from the saved .mat: one of these .mat files is 21.7 GB because it retained
% bootstrap samples, and loading it would dominate the runtime. The .nii IS the
% thresholded map - verified: PDM1.nii nonzero count equals WfullThresh.
%
% -------------------------------------------------------------------------
% by: Lukas Van Oudenhove  |  KU Leuven, September 2026
% -------------------------------------------------------------------------

analyses = {
 'discoverie model_2f, whole GM, nocov',        '/data/proj_discoverie/secondlevel/model_2f_combat_conditions/results/mediation_analysis/pdm/stressVsControl'
 'discoverie model_2f, neurosynth mask, nocov', '/data/proj_discoverie/secondlevel/model_2f_combat_conditions/results/mediation_analysis/pdm/stressVsControl_nsmask'
 'cfs model_2a, whole GM, unadjusted',          '/data/proj_cfs/secondlevel/model_2a_casecontrol_cov_scanner/results/mediation_analysis/pdm/stressVsControl'
 'cfs model_2a, whole GM, adj scanner',         '/data/proj_cfs/secondlevel/model_2a_casecontrol_cov_scanner/results/mediation_analysis/pdm/stressVsControl_adj_scanner'
 'cfs model_2a, neurosynth mask, adj scanner',  '/data/proj_cfs/secondlevel/model_2a_casecontrol_cov_scanner/results/mediation_analysis/pdm/stressVsControl_adj_scanner_nsmask'
};

% Only report analyses belonging to THIS study. rootdir is set by the model's
% s0, so running the script from each superdataset publishes each study's
% figures into its own results/html rather than dumping both studies into
% whichever project the script happened to be launched from.
if exist('rootdir','var') && ~isempty(rootdir)
    keep = startsWith(analyses(:,2), rootdir);
    fprintf('\nrootdir = %s\n  %d of %d analyses belong to this study\n', ...
        rootdir, sum(keep), numel(keep));
    analyses = analyses(keep, :);
end
if isempty(analyses)
    fprintf('\nno analyses for this study; nothing to do\n');
    return
end

fprintf('\nloading atlas for region labelling...\n');
atl = load_atlas('canlab2024_fine_2mm');

PDMREPORT = struct('analysis',{},'pdm',{},'nsig',{},'nregions',{},'table',{});

for a = 1:size(analyses,1)

    fprintf('\n\n');
    printhdr(sprintf('PDM REPORT: %s', analyses{a,1}));
    fprintf('\n\n');

    condir = analyses{a,2};
    if ~exist(condir,'dir')
        fprintf('  directory not found, skipping: %s\n', condir);
        continue
    end

    for k = 1:3
        fn = fullfile(condir, sprintf('PDM%d.nii', k));
        if ~exist(fn,'file')
            fprintf('  %s not found, skipping\n', fn);
            continue
        end
        d = fmri_data(fn, 'noverbose');
        label = sprintf('PDM%d - %s', k, analyses{a,1});
        rep = LaBGAScore_pdm_report_image(d, label, 'atlas', atl);

        PDMREPORT(end+1) = struct('analysis', analyses{a,1}, 'pdm', k, ...
            'nsig', rep.nsig, 'nregions', numel(rep.region), 'table', {rep.table}); %#ok<SAGROW>
    end

end

%% SUMMARY ACROSS ANALYSES
fprintf('\n\n');
printhdr('PDM REPORT SUMMARY');
fprintf('\n\n');
fprintf('  %-46s %5s %8s %10s\n', 'analysis', 'PDM', 'voxels', 'regions');
for a = 1:numel(PDMREPORT)
    fprintf('  %-46s %5d %8d %10d\n', PDMREPORT(a).analysis, PDMREPORT(a).pdm, ...
        PDMREPORT(a).nsig, PDMREPORT(a).nregions);
end
