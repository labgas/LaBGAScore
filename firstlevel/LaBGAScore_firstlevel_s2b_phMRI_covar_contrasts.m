%% LaBGAScore_firstlevel_s2b_phMRI_covar_contrasts.m
%
%
% *USAGE*
%
% This scripts performs first-level model contrast definition and estimation 
% for covariance analysis with a set of covariates, including spline
% interpolation
%
% It needs to be run after first-level model definition and estimation using 
% LaBGAScore_firstlevel_s1b_fit_phMRI_model.m
%
%
% *OPTIONS*
%
% For more info on the CANlab-style DSGN structure and phDSGN fields, see
% LaBGAScore_firstlevel_s1b_fit_phMRI_model.m, from which phDSGN should be
% copied (see "SET COVARIATE INFO" section below)
%
% Covariate options (STUDY-SPECIFIC, set in the "SET COVARIATE INFO"
% section below):
%
% * name_covars       cell array of covariate column names, as they appear in covars_filename
% * covars_filename   name of the .csv file with covariate values, located in BIDSdir
%
% The contrast weights are placed by column index, taken from SPM.Sess in
% each subject's SPM.mat, so the number of noise regressors per session
% (which varies when spikes are modelled) does not need to be specified.
%
% Covariate contrasts are appended to the existing ones; the script errors
% if a subject's SPM.mat already contains contrasts for these covariates.
%
% *DEPENDENCIES*
%
% SPM12 on Matlab path (spm_jobman)
% git-annex on the system path, only if SPM.mat is already annexed (the
% script then runs `git annex unlock` on SPM.mat and DSGN.mat before writing)
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
% date:   KU Leuven, August, 2024
%
% -------------------------------------------------------------------------
%
% LaBGAScore_firstlevel_s2b_phMRI_covar_contrasts.m         v1.3
%
% last modified: 2026/10/30
%
%
%% INITIATE DSGN STRUCTURE
% -------------------------------------------------------------------------

% NOTE: STUDY-SPECIFIC, to be copied from LaBGAScore_firstlevel_s1b_fit_phMRI_model.m

% DSGN INFO

phDSGN = struct();

% non-phMRI specific
phDSGN.modelingfilesdir = 'model_1_basic';
phDSGN.tr = 2.5;               % TR in seconds
phDSGN.t = 46;                 % microtime resolution, set to nr of slices is slice timing is performed - see spm documentation for more info
phDSGN.t0 = 23;                %
phDSGN.hpf = 3540;             % length of high pass filter in seconds, very long in case of phMRI
phDSGN.multireg = 'noise_regs';

% phMRI-specific
phDSGN.nr_dyns = 1416;         % nr of dynamics per run
phDSGN.start_dyn = 264;        % nr of dynamic where administration starts
phDSGN.timebin_length = 60;    % in seconds
phDSGN.sessions = {'Sucrose','Erythritol','Water'}; % names of conditions AS THEY APPEAR IN COVARIATE DATA FILE - in order of corresponding sessions!

% calculations, don't change
phDSGN.timebin_length_tr = phDSGN.timebin_length/phDSGN.tr;
phDSGN.nr_timebins = (phDSGN.nr_dyns-phDSGN.start_dyn)/phDSGN.timebin_length_tr; % nr of post-infusion timebins


%% SET COVARIATE INFO
% -------------------------------------------------------------------------

% NOTE: STUDY-SPECIFIC

name_covars = {'delta_GLP1','delta_PYY','delta_insulin','delta_leptin'};
covars_filename = 'blood_hormones_msd.csv';


%% GET & SET PATHS
% -------------------------------------------------------------------------

LaBGAScore_prep_s0_define_directories;

firstleveldir = fullfile(rootdir, 'firstlevel');

    if ~exist (firstleveldir,'dir')
        error('\nfirstlevel subdataset %s does not exist, please create using datalad commands prior to proceeding\n',firstleveldir);
    else
        firstlevelmodeldir = fullfile(firstleveldir,phDSGN.modelingfilesdir);
        if ~exist (firstlevelmodeldir,'dir')
            mkdir(firstlevelmodeldir);
        end
        phDSGN.modeldir = firstlevelmodeldir;
    end
    
phDSGN.subjects = derivsubjdirs';


%% INTERPOLATE VARIABLES
% -------------------------------------------------------------------------

num_covars = size(name_covars,2);
num_sess = size(phDSGN.sessions,2);
num_subs = size(phDSGN.subjects,2);

covars_table = readtable(fullfile(BIDSdir,covars_filename));
covars_table.Condition = categorical(covars_table.Condition);
covars_table.Timepoint(covars_table.Timepoint == -10) = 0; % recode to 0 for correct interpolation
covars_table.Scan_included = logical(covars_table.Scan_included);

missing_sessions = setdiff(phDSGN.sessions,categories(covars_table.Condition));
    if ~isempty(missing_sessions)
        error('\nsession name(s) %s not found in the Condition column of %s (capitalization matters), please check phDSGN.sessions before proceeding\n',strjoin(missing_sessions,', '),covars_filename);
    end

covar_sessions = phDSGN.sessions; % kept separately, since loading each subject's DSGN.mat below overwrites phDSGN

covars = struct();
    for s = 1:num_sess
        covars.(phDSGN.sessions{s}) = cell(num_covars,num_subs);
    end
    
subs2include = [];
idx_2include = false(1,num_subs);

for sub = 1:num_subs
    
    sub_nr = phDSGN.subjects{sub}(end-2:end);
    
    covars_table_sub = covars_table(contains(covars_table.ID,sub_nr),:);
    covars_table_sub = covars_table_sub(covars_table_sub.Scan_included,:);
    
    for c = 1:num_covars
    
        if sum(~isnan(covars_table_sub.(name_covars{c}))) < 6 % hormones only available for one session, hence no contrasts can be made; assuming covars are either all available, or all missing
            if c == 1
                idx_2include(1,sub) = 0;
            end
            continue
        else
            if c == 1
                idx_2include(1,sub) = 1;
                subs2include = [subs2include; ['sub-' sub_nr]];
            end
                
            for s = 1:num_sess
                covar_time_2interp = covars_table_sub.Timepoint(covars_table_sub.Condition == phDSGN.sessions{s});
                covar_2interp = covars_table_sub.(name_covars{c})(covars_table_sub.Condition == phDSGN.sessions{s});
                
                    if sum(~isnan(covar_2interp)) > 2
                        covar_interp_sess = interp1(covar_time_2interp,covar_2interp,[1:phDSGN.nr_timebins],'spline');
                    else
                        covar_interp_sess = [];
                    end
                    
                covars.(phDSGN.sessions{s}){c,sub} = covar_interp_sess;
            end
          
        end
        
    end
    
end 

num_subs_included = size(subs2include,1);

for s = 1:num_sess
    covars.(phDSGN.sessions{s}) = covars.(phDSGN.sessions{s})(:,idx_2include);
end


%% CREATE, SAVE, AND RUN CONTRAST BATCH
% -------------------------------------------------------------------------

% Weights are placed by column index from SPM.Sess, so they land on the
% intended timebin regressors whatever the number of noise regressors per
% session. Contrasts are named <covariate>_<cond_a>_<cond_b>_bin_<a>, for
% every pair of sessions in the subject's model (in nchoosek order), and are
% only written if the covariate is available for both sessions of the pair.

clear sub c s

for sub = 1:num_subs_included
    
    subfirstlevelmodeldir = fullfile(firstlevelmodeldir,subs2include(sub,:));
    load(fullfile(subfirstlevelmodeldir,'SPM.mat'));
    load(fullfile(subfirstlevelmodeldir,'DSGN.mat'));
    
    clear matlabbatch
    
    matlabbatch{1}.spm.stats.con.spmmat = {fullfile(subfirstlevelmodeldir,'SPM.mat')};
    
    num_existing_contrasts = size(SPM.xCon,2);
    
        if num_existing_contrasts > 0 && any(startsWith({SPM.xCon(:).name},name_covars))
            error('\nSPM.mat of %s already contains covariate contrasts for %s, delete them (or refit with s1b) before rerunning this script\n',subs2include(sub,:),strjoin(name_covars,', '));
        end
    
    num_sess_sub = size(SPM.Sess,2);
    nr_cols = size(SPM.xX.X,2);
    sess_sub_names = cell(1,num_sess_sub);
    sess_sub_covars = cell(1,num_sess_sub); % field of the covars struct matching each session
    sess_cols = cell(1,num_sess_sub);
    
    for n = 1:num_sess_sub
         sess_sub_names{n} = regexprep(char(SPM.Sess(n).U(1).name),'_bin_\d+$','');
         covar_idx = strcmpi(covar_sessions,sess_sub_names{n});
             if sum(covar_idx) ~= 1
                 error('\nsession %s of %s matches %d entries of phDSGN.sessions in this script, please check before proceeding\n',sess_sub_names{n},subs2include(sub,:),sum(covar_idx));
             end
         sess_sub_covars{n} = covar_sessions{covar_idx};
         sess_cols{n} = SPM.Sess(n).col(1:phDSGN.nr_timebins);
    end
    
        if num_sess_sub >= 2
            sess_pairs = nchoosek(1:num_sess_sub,2);
        else
            sess_pairs = zeros(0,2);
        end
    
    con = 0;
    
     for c = 1:num_covars
        
        for p = 1:size(sess_pairs,1)
            
            i = sess_pairs(p,1);
            j = sess_pairs(p,2);
            covar_i = covars.(sess_sub_covars{i}){c,sub};
            covar_j = covars.(sess_sub_covars{j}){c,sub};
            
                if isempty(covar_i) || isempty(covar_j)
                    continue
                end
            
            for a = 1:phDSGN.nr_timebins
                
                con = con + 1;
                weights = zeros(1,nr_cols);
                weights(sess_cols{i}(a)) = covar_i(1,a);
                weights(sess_cols{j}(a)) = -covar_j(1,a);
                
                matlabbatch{1}.spm.stats.con.consess{con}.tcon.name = [name_covars{c} '_' sess_sub_names{i} '_' sess_sub_names{j} '_bin_' num2str(a)];
                matlabbatch{1}.spm.stats.con.consess{con}.tcon.weights = weights;
                matlabbatch{1}.spm.stats.con.consess{con}.tcon.sessrep = 'none';
                
                phDSGN.contrastnames{num_existing_contrasts + con} = matlabbatch{1}.spm.stats.con.consess{con}.tcon.name;
                phDSGN.contrastweights{num_existing_contrasts + con} = weights;
                
            end
            
        end
        
     end
     
        if con == 0
            fprintf('\nWARNING: no covariate contrasts could be made for %s, skipping\n',subs2include(sub,:));
            continue
        end
        
    fprintf('\n%s: adding %d covariate contrasts to %d existing ones\n',subs2include(sub,:),con,num_existing_contrasts);
    
    matlabbatch{1}.spm.stats.con.delete = 0;
    
    
    % SAVE BATCHES AND RUN
    %---------------------
    
    % SPM.mat and DSGN.mat are modified in place, so unlock them first if
    % they were already annexed (datalad save turns them into read-only
    % symlinks)
    
    for f = {'SPM.mat','DSGN.mat'}
        if system(['test -L "' fullfile(subfirstlevelmodeldir,f{1}) '"']) == 0
            [status,result] = system(['cd "' subfirstlevelmodeldir '" && git annex unlock ' f{1}]);
                if status ~= 0
                    error('\ngit annex unlock %s failed in %s:\n%s\n',f{1},subfirstlevelmodeldir,result);
                end
        end
    end
    
    cd(subfirstlevelmodeldir);
    
    eval(['save cov_contrasts_first_level_' subs2include(sub,:) '.mat matlabbatch']);

    spm_jobman('initcfg');
    spm_jobman('run',matlabbatch);
    
    cd(rootdir);
    
    
    % ADD CON IMAGES NAMES TO phDSGN STRUCTURE AND SAVE
    %--------------------------------------------------
    
    % taken from SPM.xCon, so phDSGN.contrasts{c} is the image of
    % phDSGN.contrastnames{c}
    
    load(fullfile(subfirstlevelmodeldir,'SPM.mat'),'SPM');
    
    for conimg = (num_existing_contrasts+1):size(SPM.xCon,2)
        phDSGN.contrasts{conimg} = fullfile(subfirstlevelmodeldir,SPM.xCon(conimg).Vcon.fname);
    end
    
    clear SPM
    
    savefilename = fullfile(subfirstlevelmodeldir, 'DSGN.mat');
    save(savefilename, 'phDSGN', '-append');
 
    
end
