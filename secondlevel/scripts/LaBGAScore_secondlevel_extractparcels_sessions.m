%% LaBGAScore_secondlevel_extractparcels_sessions.m
%
%
% USAGE
%
% This simple script allows you to extract parcel- and roi-wise
% within-session contrasts based on the results of a between-session
% contrast and create a results table that can be flexibly used for
% plotting, stats, etc
%
%__________________________________________________________________________
%
% Author: Lukas Van Oudenhove
% date: KU Leuven, April 2026    
%
%__________________________________________________________________________
% @(#)% LaBGAScore_secondlevel_extractparcels_sessions.m          v1.0
% last modified: 2026/04/07


%% ========================================================================
% 0. USER SETTINGS — EDIT THESE
% =========================================================================

% CHOOSE ROI AND/OR PARCELS, AND SIGNIFICANCE LEVEL FOR PARCELS

do_parcel = true;
    parcel_type = 'unc'; % | 'fdr' | 'Bayes'
do_roi = true;

% INPUT DIRECTORIES

moodbugs2_prep_s0_define_directories;
moodbugs2_secondlevel_m1_s0_a_set_up_paths_always_run_first;
load(fullfile(resultsdir,'image_names_and_setup.mat'));

% SESSION INFO

nr_sess = 2; 
names_sess = {'ses-01','ses-02'};
within_session_contrast_idx = [1,2]; % indices from DAT.contrastnames
between_session_contrast_idx = 5;

% SET MANDATORY OPTIONS FROM CORRESPONDING PREP_3a_SCRIPT USED FOR
% PARCELWISE AND/OR ROI ANALYSIS

mygroupnamefield = 'contrasts'; 
results_suffix = '';

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


%% ========================================================================
% 1. EXTRACT DATA
% =========================================================================


if do_parcel
    
    % PRE-ALLOCATE
    
    betas_within_sess = cell(1,nr_sess);
    results_tables_parcel = cell(1,nr_sess);
    
    % LOAD ATLAS USED IN PREP_3a
    
    atlas = load_atlas(atlasname_glm);
    atlas_downsample = downsample_parcellation(atlas,['labels_' num2str(atlas_granularity)]);

    % LOAD RESULTS

    resultsstring = 'parcelwise_stats_and_maps_';
    analysis_type = 'parcel-wise';

    fprintf('\n\n');
    printhdr('LOADING PARCELWISE RESULTS');
    fprintf('\n\n');

    savefilenamedata = fullfile(resultsdir, [resultsstring, mygroupnamefield, '_', scaling_string, '_', results_suffix, '.mat']);

    if exist(savefilenamedata,'file')
        fprintf('\nLoading %s regression results and maps from %s\n\n', analysis_type, savefilenamedata);
        load(savefilenamedata);
    else
        fprintf('\nNo saved results file %s. Skipping this analysis.', savefilenamedata);
        fprintf('\nRun prep_3a_run_second_level_regression_and_save.m to get %s regression results first.\n', analysis_type); 
        return
    end

    % EXTRACT BETAS FOR WITHIN-SESSION CONTRASTS

    for s = 1:nr_sess

    betas_within_sess{s} = parcelwise_stats_results{1,within_session_contrast_idx(1,s)}.datmatrix;

    end
    
    % EXTRACT SIGNIFICANT PARCELS FOR BETWEEN-SESSION CONTRASTS
    
    switch parcel_type
        
        case 'unc'
            region_objs = region_objs_unc;
            regions_between_sess = region_objs{1,between_session_contrast_idx}{1,1};  
        case 'fdr'
            region_objs = region_objs_fdr;
            regions_between_sess = region_objs{1,between_session_contrast_idx}{1,1};
        case 'Bayes'
            region_objs = region_objs_Bayes;
            regions_between_sess = region_objs{1,between_session_contrast_idx}{1,1};
    end
    
    % EXTRACT PARCEL NAMES AND THEIR INDICES FROM REGION OBJECTS
    
    for r = 1:size(regions_between_sess,2)
        keyword{r} = regions_between_sess(1,r).shorttitle;
        idx(r,:) = strcmp(atlas_downsample.labels,keyword{r});
    end

    idx_sum = logical(sum(idx,1));
    
    % CREATE TABLE PER SESSION

    for s = 1:nr_sess

        betas_within_sess{s} = betas_within_sess{s}(:,idx_sum);
        results_tables_parcel{s} = array2table(betas_within_sess{s});
        results_tables_parcel{s}.PPID = DAT.BEHAVIOR.behavioral_data_table.participant_id; 
        results_tables_parcel{s}.group = DAT.BEHAVIOR.behavioral_data_table.group; 
        results_tables_parcel{s}.session = s.*ones(height(results_tables_parcel{s}),1);

            for v = 1:size(regions_between_sess,2)
                results_tables_parcel{s}.Properties.VariableNames{v} = keyword{v};
            end

    end
    
    % CONCATENATE SESSION TABLES
    
    results_table_parcel = vertcat(results_tables_parcel{:});

end


if do_roi
    
    % PRE-ALLOCATE

    results_tables_roi = cell(1,nr_sess);
    
    % LOAD RESULTS

    resultsstring = 'roi_stats_';
    analysis_type = 'roi-wise';

    fprintf('\n\n');
    printhdr('LOADING PARCELWISE RESULTS');
    fprintf('\n\n');

    savefilenamedata = fullfile(resultsdir, [resultsstring, mygroupnamefield, '_', scaling_string, '_', results_suffix, '.mat']);

    if exist(savefilenamedata,'file')
        fprintf('\nLoading %s regression results and maps from %s\n\n', analysis_type, savefilenamedata);
        load(savefilenamedata);
    else
        fprintf('\nNo saved results file %s. Skipping this analysis.', savefilenamedata);
        fprintf('\nRun prep_3a_run_second_level_regression_and_save.m to get %s regression results first.\n', analysis_type); 
        return
    end 
    
    % CREATE TABLE PER SESSION

    for s = 1:nr_sess

    results_tables_roi{s} = roi_means_table{1,within_session_contrast_idx(1,s)};    
    results_tables_roi{s}.PPID = DAT.BEHAVIOR.behavioral_data_table.participant_id;
    results_tables_roi{s}.session = s.*ones(height(results_tables_roi{s}),1);

    end
    
    % CONCATENATE SESSION TABLES

    results_table_roi = vertcat(results_tables_roi{:});

end