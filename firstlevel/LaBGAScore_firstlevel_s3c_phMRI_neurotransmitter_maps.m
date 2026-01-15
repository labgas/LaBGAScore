%% LaBGAScore_firstlevel_s3c_phMRI_neurotransmitter_maps.m
%
%
% *USAGE*
%
% This scripts calculates cosine similarity with selected PET neurotransmitter maps
% (Hansen et al, Nature Neuroscience 2022) on every first-level beta
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
% LaBGAScore_firstlevel_s3c_phMRI_neurotransmitter_maps.m         v1.0
%
% last modified: 2026/01/15
%
%
%% SET STUDY, MODEL AND SESSION INFO
% -------------------------------------------------------------------------

% NOTE: STUDY-SPECIFIC, to be copied from LaBGAScore_firstlevel_s1b_fit_phMRI_model.m

study_prefix = 'ery_ph';
modelingfilesdir = 'model_1_basic';
sessions = {'sucrose','erythritol','water'}; % names of conditions AS THEY APPEAR IN FIRST-LEVEL SPM.Vbeta FIELD - in order of corresponding sessions!


%% LOAD AND SELECT HANSEN NEUROTRANSMITTER MAPS
% -------------------------------------------------------------------------

nt_maps = load_image_set('hansen22');
nt_maps = fmri_data_st(nt_maps);        % handles metadata_table more gracefully

% NOTE: STUDY-SPECIFIC

name_nts = {'D2','MOR','CB1'};                                                                                     % as in nt_maps.metadata_table.target
results_suffix = 'nt_maps';                                                                                              % arbitrary name to append to results excel file

nt_maps_target = nt_maps.get_wh_image((contains(nt_maps.metadata_table.target,name_nts)));
nt_maps_target = nt_maps_target.get_wh_image(nt_maps_target.metadata_table.N > 10);
nt_maps_target = nt_maps_target.get_wh_image(~contains(nt_maps_target.metadata_table.tracer,'RACLOPRIDE') & ~contains(nt_maps_target.metadata_table.tracer,'FLB457'));



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
    
    nt_sim = cell(height(nt_maps_target.metadata_table),1);
    stats = cell(height(nt_maps_target.metadata_table),1);
    
        for s = 1:size(nt_sim,1)
            [stats{s}, ~ , ~ , ~ , ~] = image_similarity_plot(betas_oi_obj, 'mapset', nt_maps_target.get_wh_image(s), 'cosine_similarity', 'noplot', 'nofigure', 'notable');
            nt_sim{s} = stats{s}.r;
        end
    
    % PUT RESULTS IN TABLE
    %---------------------
    
    for n = 1:height(nt_maps_target.metadata_table)
        names{n} = [nt_maps_target.metadata_table.target{n} ' ' nt_maps_target.metadata_table.primary_reference{n}];
    end
    
    subj_table = table('Size',[size(betas_oi_fp,1) 10],'VariableNames',[{'PPID','beta_number','beta_descrip','condition','timebin'},names],'VariableTypes',{'cellstr','cellstr','cellstr','categorical','double','double','double','double','double','double'});
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
    
        for nt = 1:size(nt_sim,1)
            subj_table.(subj_table.Properties.VariableNames{varcounter}) = nt_sim{nt}';
            varcounter = varcounter + 1;
        end
        
    % ADD COVARIATES TO TABLE IF REQUESTED
    %-------------------------------------
    
        if add_covars
            
            sub_nr = phDSGN.subjects{sub}(end-2:end);
    
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
    path_result = fullfile(secondlevelmodeldir, ['NeurotransmitterResults_' results_suffix '.xlsx']);
catch
    path_result = fullfile(secondlevelmodeldir, 'Neurotransmitter.xlsx');
end

writetable(master_table, path_result);

