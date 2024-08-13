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
% For more info on the CANlab-style DSGN structure, see the scripts below
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
% LaBGAScore_firstlevel_s2b_phMRI_covar_contrasts.m         v1.0
%
% last modified: 2024/08/13
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
nr_noise_reg = 25;                                                          % 24hmp + csf, no spikes modelled in phMRI
covars_filename = 'blood_hormones_msd.csv';


%% GET & SET PATHS
% -------------------------------------------------------------------------

ery_ph_prep_s0_define_directories;

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

covars = struct();
    for s = 1:num_sess
        covars.(phDSGN.sessions{s}) = cell(num_covars,num_subs);
    end
    
subs2include = [];

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
idx_2include = logical(idx_2include);

for s = 1:num_sess
    covars.(phDSGN.sessions{s}) = covars.(phDSGN.sessions{s})(:,idx_2include);
end


%% CREATE, SAVE, AND RUN CONTRAST BATCH
% -------------------------------------------------------------------------

clear sub c s

for sub = 3:num_subs_included
    
    subfirstlevelmodeldir = fullfile(firstlevelmodeldir,subs2include(sub,:));
    load(fullfile(subfirstlevelmodeldir,'SPM.mat'));
    load(fullfile(subfirstlevelmodeldir,'DSGN.mat'));
    
    clear matlabbatch
    
    matlabbatch{1}.spm.stats.con.spmmat = {fullfile(subfirstlevelmodeldir,'SPM.mat')};
    
    num_sess_sub = size(SPM.Sess,2);
    sess_sub_indices = 1:num_sess_sub;
    num_contrasts = sum(sess_sub_indices(1,1:end-1));
    
    for n = 1:num_sess_sub
         cond = strsplit(char(SPM.Sess(n).U(1).name),'_');
         condname = cond{1};
         sess_sub_names{n} = condname;
    end
    
    num_existing_contrasts = size(SPM.xCon,2);
    
     for c = 1:num_covars
        
        for cont = 1:num_contrasts
            
            for a = 1:phDSGN.nr_timebins
                
                if cont == 1
                    
                    if ~isempty(covars.(mlreportgen.utils.capitalizeFirstChar(sess_sub_names{1})){c,sub}) && ~isempty(covars.(mlreportgen.utils.capitalizeFirstChar(sess_sub_names{2})){c,sub})
                
                        matlabbatch{1}.spm.stats.con.consess{a + (c-1)*(phDSGN.nr_timebins*3)}.tcon.name = [name_covars{c} '_' sess_sub_names{1} '_' sess_sub_names{2} '_bin_' num2str(a)];
                        matlabbatch{1}.spm.stats.con.consess{a + (c-1)*(phDSGN.nr_timebins*3)}.tcon.weights = [zeros(1, a-1) covars.(mlreportgen.utils.capitalizeFirstChar(sess_sub_names{1})){c,sub}(1,a) zeros(1,(phDSGN.nr_timebins-1)) zeros(1,nr_noise_reg) -covars.(mlreportgen.utils.capitalizeFirstChar(sess_sub_names{2})){c,sub}(1,a)]; 
                        matlabbatch{1}.spm.stats.con.consess{a + (c-1)*(phDSGN.nr_timebins*3)}.tcon.sessrep = 'none';

                        phDSGN.contrastnames{a + num_existing_contrasts + (c-1)*(phDSGN.nr_timebins*3)} = matlabbatch{1}.spm.stats.con.consess{a + (c-1)*(phDSGN.nr_timebins*3)}.tcon.name;
                        phDSGN.contrastweights{a + num_existing_contrasts + (c-1)*(phDSGN.nr_timebins*3)} = matlabbatch{1}.spm.stats.con.consess{a + (c-1)*(phDSGN.nr_timebins*3)}.tcon.weights;
                        
                    end
                
                else

                    if cont == 2
                        
                        if ~isempty(covars.(mlreportgen.utils.capitalizeFirstChar(sess_sub_names{1})){c,sub}) && ~isempty(covars.(mlreportgen.utils.capitalizeFirstChar(sess_sub_names{2})){c,sub})

                            b = a + (phDSGN.nr_timebins*(cont-1));
                        
                        else
                        
                            b = a;
                        
                        end
                        
                        if ~isempty(covars.(mlreportgen.utils.capitalizeFirstChar(sess_sub_names{1})){c,sub}) && ~isempty(covars.(mlreportgen.utils.capitalizeFirstChar(sess_sub_names{3})){c,sub})

                            matlabbatch{1}.spm.stats.con.consess{b + (c-1)*(phDSGN.nr_timebins*3)}.tcon.name = [name_covars{c} '_' sess_sub_names{1} '_' sess_sub_names{3} '_bin_' num2str(a)];
                            matlabbatch{1}.spm.stats.con.consess{b + (c-1)*(phDSGN.nr_timebins*3)}.tcon.weights = [zeros(1, a-1) covars.(mlreportgen.utils.capitalizeFirstChar(sess_sub_names{1})){c,sub}(1,a) zeros(1,((phDSGN.nr_timebins*cont-1))) zeros(1,nr_noise_reg) zeros(1,nr_noise_reg) -covars.(mlreportgen.utils.capitalizeFirstChar(sess_sub_names{3})){c,sub}(1,a)];
                            matlabbatch{1}.spm.stats.con.consess{b + (c-1)*(phDSGN.nr_timebins*3)}.tcon.sessrep = 'none';

                            phDSGN.contrastnames{b + num_existing_contrasts + (c-1)*(phDSGN.nr_timebins*3)} = matlabbatch{1}.spm.stats.con.consess{b + (c-1)*(phDSGN.nr_timebins*3)}.tcon.name;
                            phDSGN.contrastweights{b + num_existing_contrasts + (c-1)*(phDSGN.nr_timebins*3)} = matlabbatch{1}.spm.stats.con.consess{b + (c-1)*(phDSGN.nr_timebins*3)}.tcon.weights;
                            
                        end

                    elseif cont == 3
                        
                        if ~isempty(covars.(mlreportgen.utils.capitalizeFirstChar(sess_sub_names{1})){c,sub}) && ~isempty(covars.(mlreportgen.utils.capitalizeFirstChar(sess_sub_names{2})){c,sub}) && ~isempty(covars.(mlreportgen.utils.capitalizeFirstChar(sess_sub_names{3})){c,sub})
                            
                            b = a + (phDSGN.nr_timebins*(cont-1));
                        
                        else
                        
                            b = a;
                            
                        end
                        
                        if ~isempty(covars.(mlreportgen.utils.capitalizeFirstChar(sess_sub_names{2})){c,sub}) && ~isempty(covars.(mlreportgen.utils.capitalizeFirstChar(sess_sub_names{3})){c,sub})

                            matlabbatch{1}.spm.stats.con.consess{b + (c-1)*(phDSGN.nr_timebins*3)}.tcon.name = [name_covars{c} '_' sess_sub_names{2} '_' sess_sub_names{3} '_bin_' num2str(a)];
                            matlabbatch{1}.spm.stats.con.consess{b + (c-1)*(phDSGN.nr_timebins*3)}.tcon.weights = [zeros(1, phDSGN.nr_timebins) zeros(1,nr_noise_reg) zeros(1,a-1) covars.(mlreportgen.utils.capitalizeFirstChar(sess_sub_names{2})){c,sub}(1,a) zeros(1,(phDSGN.nr_timebins-1)) zeros(1,nr_noise_reg) -covars.(mlreportgen.utils.capitalizeFirstChar(sess_sub_names{3})){c,sub}(1,a)]; 
                            matlabbatch{1}.spm.stats.con.consess{b + (c-1)*(phDSGN.nr_timebins*3)}.tcon.sessrep = 'none';

                            phDSGN.contrastnames{b + num_existing_contrasts + (c-1)*(phDSGN.nr_timebins*3)} = matlabbatch{1}.spm.stats.con.consess{b + (c-1)*(phDSGN.nr_timebins*3)}.tcon.name;
                            phDSGN.contrastweights{b + num_existing_contrasts + (c-1)*(phDSGN.nr_timebins*3)} = matlabbatch{1}.spm.stats.con.consess{b + (c-1)*(phDSGN.nr_timebins*3)}.tcon.weights;
                            
                        else
                            
                            break
                            
                        end

                    else

                        error('\n More than 3 sessions not yet implemented in script, please customize before proceeding\n'); 

                    end
                end
            
            end
         
        end
        
     end
    
    matlabbatch{1}.spm.stats.con.consess = matlabbatch{1}.spm.stats.con.consess(~cellfun('isempty',matlabbatch{1}.spm.stats.con.consess)); % get rid of empty cells for missing conditions, causes spm to error out
    
    matlabbatch{1}.spm.stats.con.delete = 0;
    
    
    % SAVE BATCHES AND RUN
    %---------------------
    
    cd(subfirstlevelmodeldir);
    ! git annex unannex *.mat
    
    eval(['save cov_contrasts_first_level_' subs2include(sub,:) '.mat matlabbatch']);

    spm_jobman('initcfg');
    spm_jobman('run',matlabbatch);
    
    cd(rootdir);
    
    
    % ADD CON IMAGES NAMES TO phDSGN STRUCTURE AND SAVE
    %--------------------------------------------------
    
    tconimgs = dir(fullfile(subfirstlevelmodeldir,'con*.nii'));
    fconimgs = dir(fullfile(subfirstlevelmodeldir,'ess*.nii'));
    conimgs = [tconimgs;fconimgs];
    
    for conimg = (num_existing_contrasts+1):size(conimgs,1)
        phDSGN.contrasts{conimg} = fullfile(conimgs(conimg).folder,conimgs(conimg).name);
    end
    
    savefilename = fullfile(subfirstlevelmodeldir, 'DSGN.mat');
    save(savefilename, 'phDSGN', '-append');
 
    
end