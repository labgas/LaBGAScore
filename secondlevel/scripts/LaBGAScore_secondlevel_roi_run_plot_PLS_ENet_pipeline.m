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
% @(#)% LaBGAScore_secondlevel_roi_run_PLS_ENet_pipeline.m          v1.0
% last modified: 2026/03/07


%% ========================================================================
% 0. USER SETTINGS — EDIT THESE
% =========================================================================

% CHOOSE PIPELINE(S)

do_pls = true;
do_enet = true;


% INPUT DIRECTORIES

moodbugs2_prep_s0_define_directories;
moodbugs2_secondlevel_m1_s0_a_set_up_paths_always_run_first;


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

cons2use = 1:5; % indices from DAT.conditions or DAT.contrasts, depending on mygroupnamefield


% SET OPTIONS FOR PLS AND ENET PIPELINES

% Partial Least Squares
% help PLS_neuroimaging_pipeline for details

opts_PLS.outerK = 4;
opts_PLS.innerK = 4;
opts_PLS.nrepeats = 50;
opts_PLS.maxLV = 3;
opts_PLS.nPerm = 10000;
opts_PLS.nBoot = 10000;
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

load(fullfile(maskdir,[roi_modelname '_combined' roi_set_name]));
atlasFile = fullfile(maskdir,['combined_' roi_set_name]);
roiNames = roi_atlas.labels';


% INPUT DATA - LOAD TABLES WITH ROI DATA AND PREP FUNCTION INPUT

% pre-allocate cell arrays to store output

roi_PLSDA_results = cell(size(cons2use,2),1);
roi_ENet_results = cell(size(cons2use,2),1);

% load results

load(fullfile(resultsdir, ['roi_stats_', mygroupnamefield, '_', scaling_string, '_', results_suffix, '.mat']));

for c = cons2use

K1_ROI_file = fullfile(resultsdir,"K1_ROI.xlsx");
K1_ROI_data = readtable(K1_ROI_file);
K1_ROI_data = K1_ROI_data(:,[1:2,63:end]); % hardcoded, can be done more elegantly - ideally no need to trim if we produce a clean file from the excel results files, adding group

% convert group var into numerical if needed

groupIsCell = iscell(K1_ROI_data.Group);

if groupIsCell

    allAreChar = all(cellfun(@ischar,K1_ROI_data));

    if allAreChar
        for i = 1:height(K1_ROI_data)
            if contains(K1_ROI_data.Group(i),'CFS')
                K1_ROI_data.GroupNum(i) = 1;
            elseif contains(K1_ROI_data.Group(i),'HC')
                K1_ROI_data.GroupNum(i) = -1;
            else
                error('wrong label')
            end
        end

    end

end

% remove rows with missings if needed

K1_ROI_data = rmmissing(K1_ROI_data);

% define X and Y for model

K1_ROI_X = table2array(K1_ROI_data(:,3:end));
K1_ROI_Y = K1_ROI_data.Group;


% OUTPUT DIRECTORY

cd(resultsdir)
mkdir("K1_ROI");


%% ========================================================================
% 1. CALL PIPELINE AND PLOTTING FUNCTIONS
% =========================================================================

% PLS

if do_pls

    K1_ROI_PLSDA_results = PLSDA_neuroimaging_pipeline(K1_ROI_X,K1_ROI_Y,opts_PLS);
    
    K1_ROI_PLSDA_table = plot_PLSDA_diagnostics_neuroimaging(K1_ROI_results, [], roiNames, atlasFile, ...
        'LV',1,'TopN',14,'TopK',14,'VIP_thresh',0.8,'stab_thresh',1.5,'MapPrctile',70,'OutPrefix','K1_ROI_PLSDA','RelaxIfEmpty',true);

end

% ELASTIC NET

if do_enet

    K1_ROI_ENet_results = ENet_neuroimaging_pipeline(K1_ROI_X,K1_ROI_Y,opts_ENet);
    
    plot_ENet_diagnostics_neuroimaging(K1_ROI_ENet_results, K1_ROI_X, K1_ROI_Y, roiNames, atlasFile, ...
        'TopN',14,'TopK',14,'FreqThresh',0.5,'WeightThresh',0,'MapPrctile',70,'DoPostSelection',true,'OutPrefix','ENet_K1_ROIs','RelaxIfEmpty',true);

end

end