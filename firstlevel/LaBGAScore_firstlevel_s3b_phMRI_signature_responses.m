%% LaBGAScore_firstlevel_s3_phMRI_signature_responses.m
%
%
% *USAGE*
%
% This scripts calculates signature response on every first-level condition
% (one per timebin per condition) and writes a summary table in long format 
% for mixed model analysis in a statistical program of your choice
%
% It needs to be run after first-level model definition and estimation using 
% LaBGAScore_firstlevel_s1b_fit_phMRI_model.m
%
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
% LaBGAScore_firstlevel_s3b_phMRI_signature_responses.m         v1.0
%
% last modified: 2025/06/24
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

    if ~exist (firstleveldir,'dir')
        error('\nfirstlevel subdataset %s does not exist, please create using datalad commands prior to proceeding\n',firstleveldir);
    else
        secondlevelmodeldir = fullfile(secondleveldir,modelingfilesdir);
        if ~exist (secondlevelmodeldir,'dir')
            mkdir(secondlevelmodeldir);
        end
    end

num_sess = size(sessions,2);
num_subs = size(derivsubjdirs,1);


%% CALCULATE SIGNATURE RESPONSE FOR BETAS AND CREATE RESULTS TABLE
% -------------------------------------------------------------------------

master_table = table();

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
    
    subj_table = table('Size',[size(betas_oi_fp,1) 10],'VariableNames',[{'PPID','beta_number','beta_descrip','condition','timebin'},name_sigs],'VariableTypes',{'cellstr','cellstr','cellstr','categorical','double','double','double','double','double','double'});
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

