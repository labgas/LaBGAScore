%% LaBGAScore_secondlevel_roi_run_plot_PLS_ENet_pipeline.m
%
%
% *USAGE*
%
% This script serves as a simple wrapper to run the PLS-DA and Elastic Net
% neuroimaging pipeline functions and their plotting functions on fMRI ROI
% data saved by the prep_3a_run_second_level_regression_and_save script.
%
% To control for covariates, list them in covariate_names below. They are
% regressed out INSIDE every cross-validation fold, with the nuisance
% coefficients estimated on the training fold only. Do NOT residualize the ROI
% data beforehand: fitting the nuisance model on the full sample uses test-fold
% information and produces a biased estimate (measurably so -- on data where a
% covariate drives all the apparent signal it pushes performance far BELOW
% chance rather than above it).
%
% For more info, type the following in your Matlab command window
%
% help PLSDA_neuroimaging_pipeline
% help plot_PLSDA_diagnostics_neuroimaging
% help ENet_neuroimaging_pipeline
% help plot_ENet_diagnostics_neuroimaging
%
% or check the READMEs in the LaBGAScore Github repo
%
% https://github.com/labgas/LaBGAScore/blob/main/secondlevel/README_PLSDA_neuroimaging_pipeline.md
% https://github.com/labgas/LaBGAScore/blob/main/secondlevel/README_PLSDA_plotting.md
% https://github.com/labgas/LaBGAScore/blob/main/secondlevel/README_ENet_neuroimaging_pipeline.md
% https://github.com/labgas/LaBGAScore/blob/main/secondlevel/README_ENet_plotting.md
%
%
% *OPTIONS*
%
% * do_pls                 run the PLS-DA pipeline, default true
% * do_enet                run the Elastic Net pipeline, default false
% * group_ID                name of the group-membership variable in the saved roi_stats_*.mat table
% * covariate_names         cellstr of column names in the saved roi_stats_*.mat table to regress out
%                           fold-wise; {} for none. These columns are excluded from the feature matrix.
% * mygroupnamefield        'conditions' or 'contrasts', must match corresponding prep_3a script
% * results_suffix          suffix used when saving prep_3a results, must match corresponding prep_3a script
% * myscaling_glm           'raw' | 'scaled' | 'scaled_contrasts', must match corresponding prep_3a script (or a2_set_default_options)
% * roi_modelname           prefix of the saved roi files from the corresponding LaBGAScore_atlas_rois_from_atlas script
% * roi_set_name            descriptive name for the set of rois, from the corresponding LaBGAScore_atlas_rois_from_atlas script
% * cons2analyze            indices from DAT.conditions or DAT.contrasts (depending on mygroupnamefield) to run the pipeline(s) on
% * opts_PLS                struct of PLS-DA pipeline options (outerK, innerK, nRepeats, maxLV, nPerm, nBoot, learningSteps, seed) — see help PLSDA_neuroimaging_pipeline
% * opts_ENet               struct of Elastic Net pipeline options (outerK, innerK, nRepeats, alphaGrid, lambdaGrid, nPerm, nBoot, learningSteps, tuneRule, seed) — see help ENet_neuroimaging_pipeline
%
%
% *DEPENDENCIES*
%
% Requires prior output from prep_3a_run_second_level_regression_and_save.m
% and LaBGAScore_atlas_rois_from_atlas.m; calls LaBGAScore_smart_parallel_pool_setup.m
% and save_all_open_figures_smart.m
%
% -------------------------------------------------------------------------
%
% modified by: Lukas Van Oudenhove
%
% date:   KU Leuven, March 2026
%
% -------------------------------------------------------------------------
%
% LaBGAScore_secondlevel_roi_run_plot_PLS_ENet_pipeline.m          v1.4
%
% last modified: 2026/09/17


%% ========================================================================
% 0. USER SETTINGS — EDIT THESE
% =========================================================================

% CHOOSE PIPELINE(S)

do_pls = true;
do_enet = false;


% INPUT DIRECTORIES

% Remember where the study's own setup put the results, so the call below can be
% checked against it (see the guard immediately after).
resultsdir_before_setup = '';
if exist('resultsdir','var'), resultsdir_before_setup = resultsdir; end


LaBGAScore_prep_s0_define_directories;
a_set_up_paths_always_run_first;

% GUARD: did the path setup just move the output directory?
%
% The call above is meant to be replaced, in a study's copy, by that study's own
% s0 (e.g. mystudy_secondlevel_m2a_s0_a_set_up_paths_always_run_first). Left as
% the generic call, it RE-DERIVES resultsdir - typically from the FIRST-LEVEL
% model name - and silently overwrites whatever the study's setup had already
% set. Every result then lands in a different model's directory while the
% published report still goes to the right one, so the split is easy to miss.
if ~isempty(resultsdir_before_setup) && ~strcmp(resultsdir_before_setup, resultsdir)
    error(['\nPATH SETUP MOVED THE RESULTS DIRECTORY.\n\n' ...
           '  before: %s\n  after : %s\n\n' ...
           'The generic a_set_up_paths_always_run_first re-derived resultsdir and\n' ...
           'discarded the one your study setup had set. In your copy of this script,\n' ...
           'replace that call with your study''s own s0 path script.\n'], ...
           resultsdir_before_setup, resultsdir);
end

load(fullfile(resultsdir,'image_names_and_setup.mat'));

group_ID = 'group'; % name of variable indicating group membership in ['roi_stats_', mygroupnamefield, '_', scaling_string, '_', results_suffix, '.mat']

covariate_names = {}; % e.g. {'age','sex'}; column names in the same table to regress out fold-wise
                      % these columns are excluded from the feature matrix (see below)


% SET MANDATORY OPTIONS FROM CORRESPONDING PREP_3a_SCRIPT

mygroupnamefield = 'contrasts'; 
results_suffix = '';
myscaling_glm = 'raw'; % if not specific in corresponding prep_3a script, get from a2_set_default_options

switch myscaling_glm

    case 'raw'
        fprintf('\nContrast calculated on raw (unscaled) condition images used in second-level GLM\n\n');
        scaling_string = 'no_scaling';

    case 'scaled'
        fprintf('\nContrast calculated on z-scored condition images used in second-level GLM\n\n');
        scaling_string = 'scaling_z_score_conditions';

    case 'scaled_contrasts'
        fprintf('\nl2norm scaled contrast images used in second-level GLM\n\n');
        scaling_string = 'scaling_l2norm_contrasts';

    otherwise
        error('\nInvalid option "%s" defined in myscaling_glm variable in a2_set_default_options script, choose between "raw", "scaled", or "scaled_constrast" given option "%s" defined in mygroupnamefield variable\n\n', myscaling_glm, mygroupnamefield);

end


% GET NAMES FROM THE CORRESPONDING LaBGASCORE_ATLAS_ROIS_FROM_ATLAS SCRIPT

roi_modelname = 'model_1_basic';          % prefix which will be added to names of saved roi files which will be written in model-specific maskdir
roi_set_name = 'MIST';                    % descriptive name for set of rois which will be included in filename


% SET CONDITIONS/CONTRASTS ON WHICH YOU WANT TO RUN THE PIPELINE(S)

cons2analyze = 1:5; % indices from DAT.conditions or DAT.contrasts, depending on mygroupnamefield


% INPUT DATA - LOAD TABLES WITH ROI DATA AND PREP FUNCTION INPUT

% pre-allocate

X_vars = cell(1, size(cons2analyze,2));
p = cell(1, size(cons2analyze,2));

% load results file

load(fullfile(resultsdir, ['roi_stats_', mygroupnamefield, '_', scaling_string, '_', results_suffix, '.mat']));
input_data = roi_means_table;
% Drop the OUTCOME by name, not by position. The old code assumed the group
% variable was positionally last, which holds only when prep_3a wrote no
% covariate columns after it. In proj_cfs the table is
% [8 ROIs, group, scanner], so 1:end-1 stripped 'scanner' and left 'group'
% itself among the features - the model would then predict group from group
% and report near-perfect accuracy.
varnames = input_data{1}.Properties.VariableNames;
if ~ismember(group_ID, varnames)
    error('group_ID ''%s'' not found in the roi_stats table (vars: %s).', ...
        group_ID, strjoin(varnames, ', '));
end
varnames = setdiff(varnames, group_ID, 'stable');

% Optional explicit feature list. Without it every remaining column is treated
% as a feature, which silently includes any covariate column that prep_3a
% appended (e.g. 'scanner') unless it is also named in covariate_names.
if exist('feature_names','var') && ~isempty(feature_names)
    missing_f = setdiff(feature_names, varnames);
    if ~isempty(missing_f)
        error('feature_names not found (or is the outcome): %s', strjoin(missing_f, ', '));
    end
    varnames = feature_names;
end

% Drop the covariate columns from the FEATURE list by name. The 1:end-1 above
% only strips the group variable, which is positionally last; covariate columns
% sit among the ROI columns, so without this they would silently be modelled as
% features as well as being regressed out.
if ~isempty(covariate_names)
    missing = setdiff(covariate_names, input_data{1}.Properties.VariableNames);
    if ~isempty(missing)
        error('covariate_names not found in roi_stats table: %s', strjoin(missing, ', '));
    end
    varnames = setdiff(varnames, covariate_names, 'stable');
end

% Hard guard on the leakage bug this script was rewritten to fix: whatever
% route built varnames above, the outcome and every covariate must be absent
% from the final feature list. The original code dropped the outcome by
% POSITION (1:end-1), which is only correct when prep_3a wrote no covariate
% column after it - where it did, the model was handed the group variable as a
% predictor of itself and reported near-perfect accuracy. Asserted rather than
% assumed, and echoed to the log so the features are visible in the report.
assert(~ismember(group_ID, varnames), ...
    'LEAKAGE: outcome ''%s'' is still in the feature list.', group_ID);
if ~isempty(covariate_names)
    leaked = intersect(covariate_names, varnames);
    assert(isempty(leaked), 'LEAKAGE: covariate(s) %s still in the feature list.', strjoin(leaked, ', '));
end
fprintf('\n%d features: %s\n', numel(varnames), strjoin(varnames, ', '));
fprintf('outcome: %s (not among the features)\n', group_ID);
if isempty(covariate_names)
    fprintf('fold-wise covariates: none\n\n');
else
    fprintf('fold-wise covariates: %s\n\n', strjoin(covariate_names, ', '));
end

% create cell arrays with X vars, and define single Y var

for x = 1:size(X_vars,2)
    X_vars{x} = table2array(input_data{x}(:,varnames));
    p{x} = size(X_vars{x},2); % number of features
end

Y_var = input_data{1}.(group_ID);

% covariate matrix, one row per subject, same row order as X_vars
if isempty(covariate_names)
    covariates = [];
else
    covariates = table2array(input_data{1}(:,covariate_names));
end


% SET OPTIONS FOR PLS AND ENET PIPELINES

% Partial Least Squares
% help PLS_neuroimaging_pipeline for details

opts_PLS.outerK = 4;
opts_PLS.innerK = 4;
opts_PLS.nRepeats = 50;   % NOTE: capital R. This was 'nrepeats' before, which the
                          % pipeline never reads, so the option silently did nothing.
opts_PLS.maxLV = 3;
opts_PLS.nPerm = 5000;
opts_PLS.nBoot = 5000;
opts_PLS.learningSteps = 6;
opts_PLS.seed = 1;
opts_PLS.covariates = covariates;
opts_PLS.covariateNames = covariate_names;

% Elastic Net
% help ENet_neuroimaging_pipeline for details

opts_ENet.outerK = 4;
opts_ENet.innerK = 4;
opts_ENet.nRepeats = 50;  % NOTE: capital R, see above
opts_ENet.alphaGrid = [0.05 0.1 0.25 0.5 0.75 0.9 1];
opts_ENet.lambdaGrid = logspace(-3,1,25);
opts_ENet.nPerm = 5000;
opts_ENet.nBoot = 5000;
opts_ENet.learningSteps = 6;
opts_ENet.tuneRule = '1se';   % '1se' | 'max', see help selectENetHyperparams
opts_ENet.seed = 1;
opts_ENet.covariates = covariates;
opts_ENet.covariateNames = covariate_names;


% LOAD ATLAS

load(fullfile(maskdir,[roi_modelname '_combined' roi_set_name '.mat']));
roiatlasFile = fullfile(maskdir,['combined_' roi_set_name '.nii']);
roiNames = roi_atlas.labels';

% T1 UNDERLAY FOR PLOTTING

T1 = which('fmriprep20_template.nii');
T1_obj = fmri_data(T1);
T1_obj_resample = resample_space(T1_obj,roi_atlas);
T1_obj_resample.write('fname',fullfile(maskdir,'fmriprep20_template_downsample.nii'),'overwrite');
T1_downsample = fullfile(maskdir,'fmriprep20_template_downsample.nii');


% OUTPUT DIRECTORY

pipeline_resultsdir = fullfile(resultsdir,'pls_enet_pipeline');

    if ~exist(pipeline_resultsdir,'dir')
        mkdir(pipeline_resultsdir);
    end

    
% START PARALLEL POOL SMARTLY

LaBGAScore_smart_parallel_pool_setup;


%% ========================================================================
% 1. CALL PIPELINE AND PLOTTING FUNCTIONS
% =========================================================================

% PLS

if do_pls
    
    PLS_results = cell(1, size(cons2analyze,2));
    PLS_tables = cell(1, size(cons2analyze,2));
    
    for idx = 1:numel(cons2analyze)

        d = cons2analyze(idx); % d = actual condition/contrast index into DAT; idx = position in cons2analyze, used to index X_vars/p/PLS_results/PLS_tables (which are built by position, not by value)

        switch mygroupnamefield

            case 'conditions'

                pipeline_resultssubdir = fullfile(pipeline_resultsdir,DAT.conditions{d});

                    if ~exist(pipeline_resultssubdir,'dir')
                        mkdir(pipeline_resultssubdir);
                    end

                cd(pipeline_resultssubdir);

            case 'contrasts'

                pipeline_resultssubdir = fullfile(pipeline_resultsdir,DAT.contrastnames{d});

                    if ~exist(pipeline_resultssubdir,'dir')
                        mkdir(pipeline_resultssubdir);
                    end

                cd(pipeline_resultssubdir);

        end


        PLS_results{idx} = PLSDA_neuroimaging_pipeline(X_vars{idx},Y_var,opts_PLS);

        [max_varY, idx_LV] = max(PLS_results{idx}.varExplainedY);

        fprintf('\nPlotting latent variable %d explaining %.2f%% of the variance in Y\n\n', idx_LV, max_varY*100);

        PLS_tables{idx} = plot_PLSDA_diagnostics_neuroimaging(PLS_results{idx}, [], roiNames, roiatlasFile, ...
            'LV',idx_LV,'TopN',min(p{idx},20),'VIP_thresh',0.8,'stab_thresh',1.5,'MapPrctile',70,'OutPrefix',[num2str(d) '_PLS'],'RelaxIfEmpty',false,'UnderlayFile',T1_downsample);

        save_all_open_figures_smart(pipeline_resultssubdir,[num2str(d) '_PLS'],{'fig','svg'},true);

        clear pipeline_resultssubdir

    end
    
    saveplsfilename = fullfile(pipeline_resultsdir,'PLS_DA.mat');
    save(saveplsfilename, 'cons2analyze','varnames','PLS_results','PLS_tables','-v7.3');

end

% ELASTIC NET

if do_enet
    
    ENet_results = cell(1, size(cons2analyze,2));
    ENet_tables = cell(1, size(cons2analyze,2));
    
    for idx = 1:numel(cons2analyze)

        d = cons2analyze(idx); % d = actual condition/contrast index into DAT; idx = position in cons2analyze, used to index X_vars/p/ENet_results/ENet_tables (which are built by position, not by value)

        switch mygroupnamefield

            case 'conditions'

                pipeline_resultssubdir = fullfile(pipeline_resultsdir,DAT.conditions{d});

                    if ~exist(pipeline_resultssubdir,'dir')
                        mkdir(pipeline_resultssubdir);
                    end

                cd(pipeline_resultssubdir);

            case 'contrasts'

                pipeline_resultssubdir = fullfile(pipeline_resultsdir,DAT.contrastnames{d});

                    if ~exist(pipeline_resultssubdir,'dir')
                        mkdir(pipeline_resultssubdir);
                    end

                cd(pipeline_resultssubdir);

        end

        opts_ENet.selectionTopK = min(20, max(3, ceil(0.25 * p{idx})));

        ENet_results{idx} = ENet_neuroimaging_pipeline(X_vars{idx},Y_var,opts_ENet);

        ENet_tables{idx} = plot_ENet_diagnostics_neuroimaging(ENet_results{idx}, X_vars{idx},Y_var, roiNames, roiatlasFile, ...
            'TopN',min(p{idx},20),'FreqThresh',0.5,'WeightThresh',0,'MapPrctile',70,'DoPostSelection',true,'OutPrefix',[num2str(d) '_ENet'],'RelaxIfEmpty',false,'UnderlayFile',T1_downsample);

        save_all_open_figures_smart(pipeline_resultssubdir,[num2str(d) '_ENet'],{'fig','svg'},true);

        clear pipeline_resultssubdir

    end
    
    saveenetfilename = fullfile(pipeline_resultsdir,'ENet.mat');
    save(saveenetfilename, 'cons2analyze','varnames','ENet_results','ENet_tables','-v7.3');

end

cd(rootdir)