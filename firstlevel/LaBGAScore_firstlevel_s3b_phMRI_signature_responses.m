%% LaBGAScore_firstlevel_s3b_phMRI_signature_responses.m
%
%
% *USAGE*
%
% This scripts calculates signature response on every first-level condition
% (one per timebin per condition) and writes a summary table in long format 
% for mixed model analysis in a statistical program of your choice.
%
% It also contains an option to interpolate (behavioral, physiological)
% covariates using splines and add the interpolated values to the same
% table
%
% It needs to be run after first-level model definition and estimation using
% LaBGAScore_firstlevel_s1b_fit_phMRI_model.m
%
%
% *OPTIONS*
%
% All STUDY-SPECIFIC, set in the sections below:
%
% * study_prefix       prefix of your study-specific <study_prefix>_prep_s0_define_directories.m script, called dynamically via eval()
% * modelingfilesdir   name of subfolder within firstlevel dir where first-level model output was saved
% * sessions           cell array of condition names, in order of corresponding sessions, as they appear in the first-level SPM.Vbeta descrip field
% * path_sigs          cell array of signature .nii paths (which(...)) or load_image_set() keywords to apply
% * name_sigs          cell array of names for every individual signature in path_sigs (expand load_image_set sets into their constituent names)
% * results_suffix     arbitrary name appended to the results .xlsx filename
% * add_covars         true/false, whether to interpolate and add covariates to the results table
% * name_covars        cell array of covariate column names, as they appear in covars_filename (only used if add_covars is true)
% * covars_filename    name of the .csv file with covariate values, located in BIDSdir (only used if add_covars is true)
%
% *DEPENDENCIES*
%
% CANlabCore functions load_image_set, apply_mask, fmri_data
%
% -------------------------------------------------------------------------
%
% Based on the following scripts 
% 
% https://github.com/labgas/proj-fodmap-fmri/blob/main/brain_imaging_analysis/covariation%20analysis%20between%20brain%20and%20symptoms/interpolate_VAS_HV.m
% https://github.com/labgas/proj-fodmap-fmri/blob/main/brain_imaging_analysis/covariation%20analysis%20between%20brain%20and%20symptoms/First_level_HV_loop_with_VAS.m
% 
% by Jie Wu & Lukas Van Oudenhove
%
% https://github.com/labgas/LaBGAScore/blob/main/firstlevel/LaBGAScore_firstlevel_s1a_options_dsgn_multisess_multitask.m
%
% https://github.com/labgas/LaBGAScore/blob/main/firstlevel/LaBGAScore_firstlevel_s2a_fit_model_multisess_multitask.m
%
% adapted by: Lukas Van Oudenhove
%
% date:   KU Leuven, June, 2025
%
% -------------------------------------------------------------------------
%
% LaBGAScore_firstlevel_s3b_phMRI_signature_responses.m         v2.2
%
% last modified: 2026/09/02
%
%
%% SET STUDY, MODEL AND SESSION INFO
% -------------------------------------------------------------------------

% NOTE: STUDY-SPECIFIC, to be copied from LaBGAScore_firstlevel_s1b_fit_phMRI_model.m

study_prefix = 'ery_ph';
modelingfilesdir = 'model_1_basic';
sessions = {'sucrose','erythritol','water'}; % names of conditions AS THEY APPEAR IN FIRST-LEVEL SPM.Vbeta FIELD - in order of corresponding sessions!


%% SET SIGNATURE INFO
% -------------------------------------------------------------------------

% NOTE: STUDY-SPECIFIC

path_sigs = {which('PleasureSignature.nii'), which('Reward_Signature_bootstrapped_0.5.nii.gz'), load_image_set('ncs')};         % either 1) signature .nii file on Matlab path or 2) keyword from load_image_set
name_sigs = {'pleasure','reward','craving','craving drugs','craving food'};                                                     % if one of the signatures defined using load_image_set in path_sigs, list the names of all signatures included in the set separately here
results_suffix = 'reward_sigs';                                                                                                 % arbitrary name to append to results excel file


%% SET COVARIATE OPTION AND INFO
% -------------------------------------------------------------------------

% NOTE: STUDY-SPECIFIC

add_covars = true;
name_covars = {'delta_GLP1','delta_PYY','delta_insulin','delta_leptin'};
covars_filename = 'blood_hormones_msd.csv';


%% GET & SET PATHS
% -------------------------------------------------------------------------

eval([study_prefix '_prep_s0_define_directories']);

firstleveldir = fullfile(rootdir, 'firstlevel');

    if ~exist (firstleveldir,'dir')
        error('\nfirstlevel subdataset %s does not exist, please create using datalad commands prior to proceeding\n',firstleveldir);
    else
        firstlevelmodeldir = fullfile(firstleveldir,modelingfilesdir);
        if ~exist (firstlevelmodeldir,'dir')
            mkdir(firstlevelmodeldir);
        end
    end
    
firstlevelsubjdirs = dir(fullfile(firstlevelmodeldir,'sub-*'));
firstlevelsubjs = {firstlevelsubjdirs(:).name}';

secondleveldir = fullfile(rootdir, 'secondlevel');

    if ~exist (secondleveldir,'dir')
        error('\nsecondlevel subdataset %s does not exist, please create using datalad commands prior to proceeding\n',secondleveldir);
    else
        secondlevelmodeldir = fullfile(secondleveldir,modelingfilesdir);
        if ~exist (secondlevelmodeldir,'dir')
            mkdir(secondlevelmodeldir);
        end
    end

num_sess = size(sessions,2);
num_subs = size(firstlevelsubjs,1); % loop below indexes firstlevelsubjs, which need not cover every subject in derivdir


%% CALCULATE SIGNATURE RESPONSE FOR BETAS AND CREATE RESULTS TABLE
% -------------------------------------------------------------------------

master_table = table();

if add_covars
    
    num_covars = size(name_covars,2);

    covars_table = readtable(fullfile(BIDSdir,covars_filename));
    covars_table.Condition = lower(covars_table.Condition);    % will not always be needed but covariate file has capitalized condition names here
    covars_table.Condition = categorical(covars_table.Condition);
    covars_table.Timepoint(covars_table.Timepoint == -10) = 0; % recode to 0 for correct interpolation
    covars_table.Scan_included = logical(covars_table.Scan_included);
    
end

for sub = 1:num_subs
    
    % DO THE MATH
    %------------
    
    subfirstlevelmodeldir = fullfile(firstlevelmodeldir,firstlevelsubjs{sub});
    load(fullfile(subfirstlevelmodeldir,'SPM.mat'));
    load(fullfile(subfirstlevelmodeldir,'DSGN.mat'));
    
    betas = SPM.Vbeta;
    idx = contains({betas.descrip}, phDSGN.sessions);
    betas_oi = betas(idx);
    betas_oi_fp = cell(size(betas_oi,2),1);
        for b = 1:size(betas_oi_fp,1)
            betas_oi_fp{b} = fullfile(subfirstlevelmodeldir,betas_oi(b).fname(1,:));
        end
    betas_oi_obj = fmri_data(betas_oi_fp);
    
    sigs = cell(size(path_sigs,2),1);
    
        for s = 1:size(path_sigs,2)
            sigs{s} = apply_mask(betas_oi_obj,path_sigs{s},'pattern_expression');
        end
    
    % PUT RESULTS IN TABLE
    %---------------------
    
    % the first five columns are fixed; the remaining ones are one per entry in
    % name_sigs, whose number is study-specific, hence derive the table width and
    % variable types from it rather than hardcoding them
    varnames_table = [{'PPID','beta_number','beta_descrip','condition','timebin'},name_sigs];
    vartypes_table = [{'cellstr','cellstr','cellstr','categorical','double'},repmat({'double'},1,numel(name_sigs))];
    subj_table = table('Size',[size(betas_oi_fp,1) numel(varnames_table)],'VariableNames',varnames_table,'VariableTypes',vartypes_table);
    subj_table.(subj_table.Properties.VariableNames{1}) = repmat(firstlevelsubjs{sub},height(subj_table),1);
    subj_table.(subj_table.Properties.VariableNames{2}) = {betas_oi(:).fname}';
    subj_table.(subj_table.Properties.VariableNames{3}) = {betas_oi(:).descrip}';
        for d = 1:height(subj_table)
            for c = 1:size(phDSGN.sessions,2)
                if contains(subj_table.(subj_table.Properties.VariableNames{3}){d},phDSGN.sessions{c})
                   subj_table.(subj_table.Properties.VariableNames{4})(d) = phDSGN.sessions{c};
                end
            end
        end
        
    subj_table.(subj_table.Properties.VariableNames{5}) = repmat([1:phDSGN.nr_timebins]',(height(subj_table)/phDSGN.nr_timebins),1);
        
        varcounter = 6;
    
        for sig = 1:size(sigs,1)
            if size(sigs{sig},2) == 1
                subj_table.(subj_table.Properties.VariableNames{varcounter}) = sigs{sig};
                varcounter = varcounter + 1;
            else
                for subsig = 1:size(sigs{sig},2)
                    subj_table.(subj_table.Properties.VariableNames{varcounter}) = sigs{sig}(:,subsig);
                    varcounter = varcounter + 1;
                end
            end
        end
        
    % ADD COVARIATES TO TABLE IF REQUESTED
    %-------------------------------------
    
        if add_covars
            
            sub_nr = firstlevelsubjs{sub}(end-2:end); % take the subject actually being processed, not an index into phDSGN.subjects loaded from this subject's DSGN.mat
    
            covars_table_sub = covars_table(contains(covars_table.ID,sub_nr),:);
            covars_table_sub = covars_table_sub(covars_table_sub.Scan_included,:);
            
            covars = cell(num_covars,1);

            for covar = 1:num_covars
                
                covar_interp_subj = [];

                for sess = 1:(height(subj_table)/phDSGN.nr_timebins)
                    
                    covar_time_2interp = covars_table_sub.Timepoint(covars_table_sub.Condition == phDSGN.sessions{sess});
                    covar_2interp = covars_table_sub.(name_covars{covar})(covars_table_sub.Condition == phDSGN.sessions{sess});

                        if sum(~isnan(covar_2interp)) > 2
                            covar_interp_sess = interp1(covar_time_2interp,covar_2interp,[1:phDSGN.nr_timebins],'spline');
                            covar_interp_sess = covar_interp_sess';
                        else
                            covar_interp_sess = NaN(phDSGN.nr_timebins,1);
                        end

                    covar_interp_subj = [covar_interp_subj;covar_interp_sess];
                    
                end
                
                subj_table.(name_covars{covar}) = covar_interp_subj;

            end

        end
    
    master_table = [master_table;subj_table];
        
end


%% SAVE MASTER RESULTS TABLE
% -------------------------------------------------------------------------
    
try
    path_result = fullfile(secondlevelmodeldir, ['SignatureResults_' results_suffix '.xlsx']);
catch
    path_result = fullfile(secondlevelmodeldir, 'SignatureResults.xlsx');
end

writetable(master_table, path_result);

