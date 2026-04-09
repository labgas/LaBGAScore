%% LaBGAScore_secondlevel_roi_run_plot_PLS_ENet_pipeline.m
%
%
% USAGE
%
% This script serves as a simple wrapper to run the PLS-DA and Elastic Net
% neuroimaging pipeline functions and their plotting functions on fMRI ROI
% data saved by the prep_3a_run_second_level_regression_and_save script.
%
% If you want to control for covariates, residualize PRIOR to running this
% script!
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
%__________________________________________________________________________
%
% Author: Lukas Van Oudenhove
% date: KU Leuven, March 2026    
%
%__________________________________________________________________________
% @(#)% LaBGAScore_secondlevel_roi_run_PLS_ENet_pipeline.m          v1.1
% last modified: 2026/04/08


%% ========================================================================
% 0. USER SETTINGS — EDIT THESE
% =========================================================================

% CHOOSE PIPELINE(S)

do_pls = true;
do_enet = false;


% INPUT DIRECTORIES

LaBGAScore_prep_s0_define_directories;
a_set_up_paths_always_run_first;
load(fullfile(resultsdir,'image_names_and_setup.mat'));

group_ID = 'group'; % name of variable indicating group membership in ['roi_stats_', mygroupnamefield, '_', scaling_string, '_', results_suffix, '.mat']


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
varnames = input_data{1}.Properties.VariableNames(1:end-1); % group var always last

% create cell arrays with X vars, and define single Y var

for x = 1:size(X_vars,2)
    X_vars{x} = table2array(input_data{x}(:,varnames)); % group var always last
    p{x} = size(X_vars{x},2); % number of features
end

Y_var = input_data{1}.(group_ID);


% SET OPTIONS FOR PLS AND ENET PIPELINES

% Partial Least Squares
% help PLS_neuroimaging_pipeline for details

opts_PLS.outerK = 4;
opts_PLS.innerK = 4;
opts_PLS.nrepeats = 50;
opts_PLS.maxLV = 3;
opts_PLS.nPerm = 5000;
opts_PLS.nBoot = 5000;
opts_PLS.learningSteps = 6;

% Elastic Net
% help ENet_neuroimaging_pipeline for details

opts_ENet.outerK = 4;
opts_ENet.innerK = 4;
opts_ENet.nrepeats = 50;
opts_ENet.alphaGrid = [0.05 0.1 0.25 0.5 0.75 0.9 1];
opts_ENet.lambdaGrid = logspace(-3,1,25);
opts_ENet.nPerm = 5000;
opts_ENet.nBoot = 5000;
opts_ENet.learningSteps = 6;


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
    
    for d = cons2analyze
        
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
                

        PLS_results{d} = PLSDA_neuroimaging_pipeline(X_vars{d},Y_var,opts_PLS);
        
        [max_varY, idx_LV] = max(PLS_results{d}.varExplainedY); 
        
        fprintf('\nPlotting latent variable %d explaining %.2f%% of the variance in Y\n\n', idx_LV, max_varY*100);
    
        PLS_tables{d} = plot_PLSDA_diagnostics_neuroimaging(PLS_results{d}, [], roiNames, roiatlasFile, ...
            'LV',idx_LV,'TopN',min(p{d},20),'VIP_thresh',0.8,'stab_thresh',1.5,'MapPrctile',70,'OutPrefix',[num2str(cons2analyze(1,d)) '_PLS'],'RelaxIfEmpty',false,'UnderlayFile',T1_downsample);

        save_all_open_figures_smart(pipeline_resultssubdir,[num2str(cons2analyze(1,d)) '_PLS'],{'fig','svg'},true);
        
        clear pipeline_resultssubdir
        
    end
    
    saveplsfilename = fullfile(pipeline_resultsdir,'PLS_DA.mat');
    save(saveplsfilename, 'cons2analyze','varnames','PLS_results','PLS_tables','-v7.3');

end

% ELASTIC NET

if do_enet
    
    ENet_results = cell(1, size(cons2analyze,2));
    ENet_tables = cell(1, size(cons2analyze,2));
    
    for d = cons2analyze
        
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
        
        opts_ENet.selectionTopK = min(20, max(3, ceil(0.25 * p{d})));

        ENet_results{d} = ENet_neuroimaging_pipeline(X_vars{d},Y_var,opts_ENet);

        ENet_tables{d} = plot_ENet_diagnostics_neuroimaging(ENet_results{d}, X_vars{d},Y_var, roiNames, roiatlasFile, ...
            'TopN',min(p{d},20),'FreqThresh',0.5,'WeightThresh',0,'MapPrctile',70,'DoPostSelection',true,'OutPrefix',[num2str(cons2analyze(1,d)) '_ENet'],'RelaxIfEmpty',false,'UnderlayFile',T1_downsample);

        save_all_open_figures_smart(pipeline_resultssubdir,[num2str(cons2analyze(1,d)) '_ENet'],{'fig','svg'},true);

        clear pipeline_resultssubdir
    
    end
    
    saveenetfilename = fullfile(pipeline_resultsdir,'ENet.mat');
    save(saveenetfilename, 'cons2analyze','varnames','ENet_results','ENet_tables','-v7.3');

end

cd(rootdir)