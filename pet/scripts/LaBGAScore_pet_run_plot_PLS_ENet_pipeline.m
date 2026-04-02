%% LaBGAScore_pet_run_plot_PLS_ENet_pipeline.m
%
%
% USAGE
%
% This script serves as a simple wrapper to run the PLS-DA and Elastic Net
% neuroimaging pipeline functions and their plotting functions on parcel-wise 
% (whole-brain or ROI) PET data.
%
% The script takes excel outputs from the LaBGAScore_pet_model_TSPO_DPA714.m script
% as input, but residualizing for genotype (and if needed other covariates)
% needs to be done PRIOR to using this script!
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
% @(#)% LaBGAScore_pet_run_PLS_ENet_pipeline.m          v1.0
% last modified: 2026/03/07


%% ========================================================================
% 0. USER SETTINGS — EDIT THESE
% =========================================================================

% CHOOSE PIPELINE(S)

do_pls = true;
do_enet = true;

% INPUT DIRECTORIES

cfs_pet_s3_a_set_up_paths_always_run_first;


% INPUT DATA - READ FILES WITH RESIDUALIZED ROI DATA AND PREP INPUT

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

% Example 1: (selected) parcels from canlab2024

atlas = load_atlas('canlab2024');
atlas = downsample_parcellation(atlas,'labels_2');
atlas.probability_maps = [];
idx_excluded_regions = [180 181 182 183 190 191 193 201 203 209 211 214 215 216 217 218 219 220 221 222 223 229 232 245 246];
atlas_reduced = atlas.remove_atlas_region(idx_excluded_regions);
atlas_reduced.fullpath = 'canlab2024_221parcels.nii';
atlas_reduced.write
atlasFile = 'canlab2024_221parcels.nii';
roiNames = K1_parcels';

% Example 2: ROI atlas made with LaBGAScore_atlas_rois_from_atlas script

load(fullfile(maskdir,"pet_trc-DPA714_combinedTSPO.mat"));
atlasFile = 'combined_TSPO.nii';
roiNames = roi_atlas.labels';

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