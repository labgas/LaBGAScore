%% ery_ph_firstlevel_s1_fit_phMRI_model.m
%
%
% *USAGE*
%
% This scripts performs first-level model specification, estimation, and
% contrast definition and estimation for phMRI designs with two or three
% conditions (more than three would need hard-coded adaptations).
%
% *OPTIONS*
% 
% For more info on the CANlab-style DSGN and LaBGAS_options structure, see
% the scripts below
%
% -------------------------------------------------------------------------
%
% Based on the following scripts 
% 
% https://github.com/labgas/proj-fodmap-fmri/blob/main/brain_imaging_analysis/1st%20level%20and%202nd%20level%20analysis/first_level_looped_fodmap_HC.m
% by Jie Wu & Lukas Van Oudenhove
%
% https://github.com/labgas/LaBGAScore/blob/main/firstlevel/LaBGAScore_firstlevel_s1a_options_dsgn_multisess_multitask.m
%
% https://github.com/labgas/LaBGAScore/blob/main/firstlevel/LaBGAScore_firstlevel_s2a_fit_model_multisess_multitask.m
%
% adapted by: Lukas Van Oudenhove
%
% date:   KU Leuven, April, 2024
%
% -------------------------------------------------------------------------
%
% LaBGAScore_firstlevel_s1b_fit_phMRI_model.m         v1.0
%
% last modified: 2024/04/26
%
%
%% INITIATE DSGN AND LABGAS_OPTIONS STRUCTURES
% -------------------------------------------------------------------------

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
phDSGN.sessions = {'sucrose','erythritol','water'}; % names of conditions - in order of corresponding sessions!

% calculations, don't change
phDSGN.timebin_length_tr = phDSGN.timebin_length/phDSGN.tr;
phDSGN.nr_timebins = (phDSGN.nr_dyns-phDSGN.start_dyn)/phDSGN.timebin_length_tr; % nr of post-infusion timebins


% LABGAS OPTIONS

LaBGAS_options = struct();

% required options
LaBGAS_options.mandatory.omit_spike_trials = 'no';
LaBGAS_options.mandatory.spikes_percent_threshold=0.15;
LaBGAS_options.mandatory.vif_thresh=2;
LaBGAS_options.movement_reg_quadratic = true; % change to false if you don't want to add quadratic terms for movement parameters and their first-order derivatives

% spike options
LaBGAS_options.spikes.dvars_threshold = 2;
LaBGAS_options.spikes.spike_additional_vols=0; % OPTIONAL, NOT RECOMMENDED TO TURN ON


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
        phDSGN.modeldir = firstlevelmodeldir; % make CANlab DSGN structure for the purpose of using it in second level analysis
    end
    
phDSGN.subjects = derivsubjdirs';

    
%% FIT PHMRI MODEL LOOPED OVER SUBJECTS
% -------------------------------------------------------------------------

spm('defaults','FMRI')

empty_cells = @(x)isempty(x);

for sub = 1:size(derivsubjs,1)
    
    
    clear matlabbatch
    matlabbatch = struct([]); % create empty structures
    
    subfirstlevelmodeldir = fullfile(firstlevelmodeldir,derivsubjs{sub});
        if ~exist(subfirstlevelmodeldir,'dir')
            mkdir(subfirstlevelmodeldir);
        end
    
    % IMAGES AND NOISE REGRESSORS
    %--------------------------------
    
    derivimages = cell(1,size(phDSGN.sessions,2));
    fmriprep_noisefiles = cell(1,size(phDSGN.sessions,2));
    derivsubjsesdirs = cell(1,size(phDSGN.sessions,2));
    derivsubjsesmodeldirs = cell(1,size(phDSGN.sessions,2));
    noise_regs = cell(1,size(phDSGN.sessions,2));
    phDSGN.conditions = cell(1,size(phDSGN.sessions,2));

        for sess = 1:size(phDSGN.sessions,2)
            
            derivsubjsesdirs{sess} = fullfile(derivsubjdirs{sub},['ses-' num2str(sess)],'func');
            derivsubjsesmodeldirs{sess} = fullfile(derivsubjsesdirs{sess},phDSGN.modelingfilesdir);
                if exist(derivsubjsesdirs{sess},'dir')
                    if ~exist(derivsubjsesmodeldirs{sess},'dir')
                        mkdir(derivsubjsesmodeldirs{sess});
                    end
                end
                
            derivimg = dir(fullfile(derivsubjsesdirs{sess},'s6*.nii.gz'));
                if size(derivimg,1) > 1
                    error('\nMore than one image file found with filter s6*.nii.gz in dir %s, please check before proceeding\n', derivsubjsesdirs{sess});
                elseif size(derivimg,1) < 1
                    fprintf('\nWARNING: No image files found with filter s6*.nii.gz in dir %s, missing condition - proceeding\n', derivsubjsesdirs{sess});
                    continue
                else
                    derivimg_nii = fullfile(derivimg.folder, derivimg.name(1,1:end-3));
                        if ~exist(derivimg_nii,'file')
                            gunzip(fullfile(derivimg.folder,derivimg.name));
                        end
                    derivimages{sess} = spm_select('ExtFPlist', derivimg.folder, derivimg.name(1,1:end-3), Inf);
                end

            fmriprep_noisefile = dir(fullfile(derivsubjsesdirs{sess},'*desc-confounds_timeseries.tsv'));
                if size(fmriprep_noisefile,1) > 1
                    error('\nMore than one noise file found with filter *desc-confounds_timeseries.tsv in dir %s, please check before proceeding\n', derivsubjsesdirs{sess});
                elseif size(fmriprep_noisefile,1) < 1
                    fprint('\nWARNING: No noise files found with filter *desc-confounds_timeseries.tsv in dir %s, missing condition - proceeding\n', derivsubjsesdirs{sess});
                    continue
                else
                    fmriprep_noisefiles{sess} = fullfile(fmriprep_noisefile.folder,fmriprep_noisefile.name);
                end
                
        end
        
        clear sess
     
        empty_derivimgs = cellfun(empty_cells,derivimages);
        empty_noisefiles = cellfun(empty_cells,fmriprep_noisefiles);
            if ~isequal(empty_derivimgs,empty_noisefiles)
                error('\nNumber of image files does not match number of noisefiles for subject %s, please check before proceeding\n', derivsubjs{sub});
            end
            
        for sess = 1:size(phDSGN.sessions,2)
            
            if ~isempty(fmriprep_noisefiles{sess})

                % load confound regressor file generated by fMRIprep into Matlab table
                % variable
                Rfull = readtable(fmriprep_noisefiles{sess},'TreatAsEmpty','n/a','FileType', 'text', 'Delimiter', 'tab');

                % replace NaNs in first row with Os
                wh_replace = ismissing(Rfull(1,:));
                    if any(wh_replace)
                        Rfull{1, wh_replace} = zeros(1, sum(wh_replace)); % make array of zeros of the right size
                    end

                % calculate and extract confound regressors

                % define regressors in fMRIprep output
                regs=Rfull.Properties.VariableNames;
                spike_cols = contains(regs,'motion_outlier');
                Rspikes=Rfull(:,spike_cols);
                Rspikes.spikes=sum(Rspikes{:,1:end},2);
                volume_idx = [1:height(Rfull)]; 
                spikes = volume_idx(Rspikes.spikes==1);

                % flag user-specified number of volumes after each spike
                % Motion can create artifacts lasting longer than the single image we
                % usually account for using spike id scripts. we're also going to flag the
                % following TRs, the number of which is defined by the user. If
                % 'LaBGAS_options.spikes.spike_additional_vols' remains unspecified, everything will proceed as
                % it did before, meaning spikes will be identified and flagged in the
                % creation of nuisance regressors without considering the following TRs
                % Add them if user requested, for both nuisance_covs and dvars_spikes_regs
                    if LaBGAS_options.spikes.spike_additional_vols ~= 0
                        additional_spikes_regs = zeros(height(Rfull),size(spikes,2)*LaBGAS_options.spikes.spike_additional_vols);
                            % This loop will create a separate column with ones in each row (TR) 
                            % we would like to consider a nuisance regressor
                            for spike = 1:size(spikes,2) 
                                additional_spikes_regs(spikes(spike)+1 : spikes(spike)+LaBGAS_options.spikes.spike_additional_vols,(spike*LaBGAS_options.spikes.spike_additional_vols-(LaBGAS_options.spikes.spike_additional_vols-1)):(spike*LaBGAS_options.spikes.spike_additional_vols)) = eye(LaBGAS_options.spikes.spike_additional_vols);
                            end
                            clear spike
                        % if any spikes went beyond the end, trim it down
                        additional_spikes_regs = additional_spikes_regs(1:height(Rfull),:);
                        % add the additional spikes to the larger matrix
                        Rfull = [Rfull array2table(additional_spikes_regs)];
                    end

                % remove redundant spike regressors
                regs = Rfull.Properties.VariableNames;
                spike_cols = contains(regs,'motion_outlier');
                additional_spike_cols = contains(regs,'additional_spikes'); 
                [duplicate_rows, ~] = find(sum(Rfull{:, spike_cols | additional_spike_cols}, 2)>1);
                    for row = 1:length(duplicate_rows) % This loop sets duplicate values to zero; drops them later (to keep indices the same during the loop)
                        [~,curr_cols] = find(Rfull{duplicate_rows(row),:}==1);
                        Rfull{duplicate_rows(row), curr_cols(2:end)} = 0;
                    end
                    clear row
                Rfull = Rfull(1:height(Rfull), any(table2array(Rfull)));

                % Select confound and spike regressors to return for use in GLM 
                regsfull = Rfull.Properties.VariableNames;
                motion_cols = contains(regsfull,'rot') | contains(regsfull,'trans');
                motion_cols_no_quad = (contains(regsfull,'rot') | contains(regsfull,'trans')) & ~contains(regsfull,'power2');
                spike_cols = contains(regsfull,'mahal_spikes') | contains(regsfull,'motion_outlier'); 
                dvars_cols = contains(regsfull,'dvars_spikes'); 
                additional_spike_cols = contains(regsfull,'additional_spikes'); 
                    if LaBGAS_options.movement_reg_quadratic
                        Rmotion = Rfull(:,motion_cols);
                    else
                        Rmotion = Rfull(:,motion_cols_no_quad);
                    end
                Rspikes = Rfull(:,spike_cols | dvars_cols | additional_spike_cols);
                Rcsf = table(Rfull.csf,'VariableNames',{'csf'});

                % recalculate derivatives, z-score, and recalculate
                % quadratics to orthogonalize
                regsmotion = Rmotion.Properties.VariableNames;
                regsnot2keep = contains(regsmotion,'derivative') | contains(regsmotion,'power');
                regsderiv = regsmotion(~contains(regsmotion,'power'));
                regsderiv = regsderiv(contains(regsderiv,'derivative'));
                Rmotion = Rmotion(:,~regsnot2keep);
                deriv = @(x) gradient(x);
                Rmotionderiv = varfun(deriv,Rmotion);
                    for regderiv = 1:size(Rmotionderiv.Properties.VariableNames,2)
                        Rmotionderiv.Properties.VariableNames{regderiv} = regsderiv{regderiv};
                    end
                Rmotionderiv = [Rmotion,Rmotionderiv];
                regsmotion2 = Rmotionderiv.Properties.VariableNames;
                zscore = @(x) zscore(x);
                Rmotionzscore = varfun(deriv,Rmotionderiv);
                    for regz = 1:size(Rmotionzscore.Properties.VariableNames,2)
                        Rmotionzscore.Properties.VariableNames{regz} = regsmotion2{regz};
                    end

                    if LaBGAS_options.movement_reg_quadratic
                        quad = @(x) x.^ 2;
                        Rmotionquad = varfun(quad,Rmotionzscore);
                        Rmotionfinal = [Rmotionzscore,Rmotionquad];
                    else
                        Rmotionfinal = Rmotionzscore;
                    end
                R = [Rmotionfinal,Rspikes,Rcsf];

                names = R.Properties.VariableNames;

                % get row indices for spikes for later use
                Rspikes.spikes=sum(Rspikes{:,1:end},2);
                volume_idx = [1:height(Rfull)]; 
                spikes = volume_idx(Rspikes.spikes==1)';

                % compute and output how many spikes total
                n_spike_regs = sum(dvars_cols | spike_cols | additional_spike_cols);
                n_spike_regs_percent = n_spike_regs / height(Rfull);

                % print warning if #volumes identified as spikes exceeds
                % user-defined threshold
                    if n_spike_regs_percent > LaBGAS_options.mandatory.spikes_percent_threshold
                        warning('\nnumber of volumes identified as spikes exceeds threshold %s in %s',LaBGAS_options.mandatory.spikes_percent_threshold,subjrunnames{run})
                    end

                % save confound regressors as matrix named R for use in
                % SPM/CANlab GLM model tools
                R=table2array(R);
                noise_regs{sess} = R;

                % write confound regressors
                filename_noise_regs = fullfile(derivsubjsesmodeldirs{sess},'noise_regs');
                save(filename_noise_regs,'R','names');

                clear Rmotion* Rcsf Rspikes
            
            else
                continue
                
            end

        end
                

    % FIRST BATCH: MODEL SPECIFICATION
    %--------------------------------

    matlabbatch{1}.spm.stats.fmri_spec.dir = {subfirstlevelmodeldir};
    matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'scans';
    matlabbatch{1}.spm.stats.fmri_spec.timing.RT = phDSGN.tr;
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = phDSGN.t;
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = phDSGN.t0;
    
    sess_subj = phDSGN.sessions(~empty_derivimgs);
    derivimages = derivimages(~empty_derivimgs);
    derivsubjsesmodeldirs = derivsubjsesmodeldirs(~empty_derivimgs);
    
        for sess = 1:size(sess_subj,2)
            
            phDSGN.conditions{sess} = cell(1,size(sess_subj,2));
            
            matlabbatch{1}.spm.stats.fmri_spec.sess(sess).scans = cellstr(derivimages{sess});
            
            for b = 1:phDSGN.nr_timebins
                
                matlabbatch{1}.spm.stats.fmri_spec.sess(sess).cond(b).name = [sess_subj{sess} '_bin_' num2str(b)];
                matlabbatch{1}.spm.stats.fmri_spec.sess(sess).cond(b).onset = phDSGN.timebin_length_tr * b + (phDSGN.start_dyn - phDSGN.timebin_length_tr);
                matlabbatch{1}.spm.stats.fmri_spec.sess(sess).cond(b).duration = phDSGN.timebin_length_tr;
                matlabbatch{1}.spm.stats.fmri_spec.sess(sess).cond(b).tmod = 0;
                matlabbatch{1}.spm.stats.fmri_spec.sess(sess).cond(b).pmod = struct('name', {}, 'param', {}, 'poly', {});
                matlabbatch{1}.spm.stats.fmri_spec.sess(sess).cond(b).orth = 1;
                matlabbatch{1}.spm.stats.fmri_spec.sess(sess).hpf = phDSGN.hpf;
                matlabbatch{1}.spm.stats.fmri_spec.sess(sess).multi_reg = {fullfile(derivsubjsesmodeldirs{sess},[phDSGN.multireg '.mat'])};
                
                phDSGN.conditions{sess}{b} = [sess_subj{sess} '_bin_' num2str(b)];
                
            end
            
        end
    
    matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
    matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
    matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
    matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';
    matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
    matlabbatch{1}.spm.stats.fmri_spec.mthresh = -inf;

    
    % SECOND BATCH: ESTIMATION
    %------------------------
    
    matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
    matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

    
    % THIRD BATCH: CONTRASTS
    %----------------------
    
    matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    
    sess_indices = 1:size(sess_subj,2);
    n_contrasts = sum(sess_indices(1,1:end-1));
    
    if n_contrasts < 1
        fprintf('\nWARNING: Only 1 session found for subject %s, skipping contrasts\n', derivsubjs{sub});
        continue
    end
    
    % T-CONTRASTS
    
    for cont = 1:n_contrasts % coded for max 3 sessions (and hence 3 contrasts)
        
        for a = 1:phDSGN.nr_timebins %loop for contrast manager
            
            if cont == 1
                
                matlabbatch{3}.spm.stats.con.consess{a}.tcon.name = [sess_subj{1} '_' sess_subj{2} '_bin_' num2str(a)];
                matlabbatch{3}.spm.stats.con.consess{a}.tcon.weights = [zeros(1, a-1) 1 zeros(1,(phDSGN.nr_timebins-1)) zeros(1,size(noise_regs{1},2)) -1]; 
                matlabbatch{3}.spm.stats.con.consess{a}.tcon.sessrep = 'none';
                
                phDSGN.contrastnames{a} = matlabbatch{3}.spm.stats.con.consess{a}.tcon.name;
                phDSGN.contrastweights{a} = matlabbatch{3}.spm.stats.con.consess{a}.tcon.weights;
                
            else
                
                b = a + (phDSGN.nr_timebins*(cont-1));
                
                if cont == 2
                    
                    matlabbatch{3}.spm.stats.con.consess{b}.tcon.name = [sess_subj{1} '_' sess_subj{3} '_bin_' num2str(a)];
                    matlabbatch{3}.spm.stats.con.consess{b}.tcon.weights = [zeros(1, a-1) 1 zeros(1,((phDSGN.nr_timebins*cont-1))) zeros(1,size(noise_regs{1},2)) zeros(1,size(noise_regs{2},2)) -1];
                    matlabbatch{3}.spm.stats.con.consess{b}.tcon.sessrep = 'none';
                    
                    phDSGN.contrastnames{b} = matlabbatch{3}.spm.stats.con.consess{b}.tcon.name;
                    phDSGN.contrastweights{b} = matlabbatch{3}.spm.stats.con.consess{b}.tcon.weights;
                    
                elseif cont == 3
                    
                    matlabbatch{3}.spm.stats.con.consess{b}.tcon.name = [sess_subj{2} '_' sess_subj{3} '_bin_' num2str(a)];
                    matlabbatch{3}.spm.stats.con.consess{b}.tcon.weights = [zeros(1, phDSGN.nr_timebins) zeros(1,size(noise_regs{1},2)) zeros(1,a-1) 1 zeros(1,(phDSGN.nr_timebins-1)) zeros(1,size(noise_regs{3},2)) -1]; 
                    matlabbatch{3}.spm.stats.con.consess{b}.tcon.sessrep = 'none';
                    
                    phDSGN.contrastnames{b} = matlabbatch{3}.spm.stats.con.consess{b}.tcon.name;
                    phDSGN.contrastweights{b} = matlabbatch{3}.spm.stats.con.consess{b}.tcon.weights;
                    
                else
                    
                    error('\n More than 3 sessions not yet implemented in script, please customize before proceeding\n'); 
                    
                end
            end
            
        end
         
    end
    
    % F-CONTRASTS
    
    c = b+1;
    
    if n_contrasts == 1
    
        matlabbatch{3}.spm.stats.con.consess{c}.fcon.name = [sess_subj{1} '_main_effect'];
        matlabbatch{3}.spm.stats.con.consess{c}.fcon.weights = eye(phDSGN.nr_timebins);
        matlabbatch{3}.spm.stats.con.consess{c}.fcon.sessrep = 'none';
        c = c+1;
        matlabbatch{3}.spm.stats.con.consess{c}.fcon.name = [sess_subj{1} '_interaction_effect'];
        matlabbatch{3}.spm.stats.con.consess{c}.fcon.weights = diff(eye(phDSGN.nr_timebins));
        matlabbatch{3}.spm.stats.con.consess{c}.fcon.sessrep = 'none';
        c = c+1;
        matlabbatch{3}.spm.stats.con.consess{c}.fcon.name = [sess_subj{2} '_main_effect'];
        matlabbatch{3}.spm.stats.con.consess{c}.fcon.weights = [zeros(phDSGN.nr_timebins) zeros(phDSGN.nr_timebins,size(noise_regs{1},2)) eye(phDSGN.nr_timebins)];
        matlabbatch{3}.spm.stats.con.consess{c}.fcon.sessrep = 'none';
        c = c+1;
        matlabbatch{3}.spm.stats.con.consess{c}.fcon.name = [sess_subj{2} '_interaction_effect'];
        matlabbatch{3}.spm.stats.con.consess{c}.fcon.weights = [zeros(phDSGN.nr_timebins-1,phDSGN.nr_timebins) zeros(phDSGN.nr_timebins-1,size(noise_regs{1},2)) diff(eye(phDSGN.nr_timebins))];
        matlabbatch{3}.spm.stats.con.consess{c}.fcon.sessrep = 'none';
        c = c+1;
        matlabbatch{3}.spm.stats.con.consess{c}.fcon.name = [sess_subj{1} '_' sess_subj{2} '_main_effect'];
        matlabbatch{3}.spm.stats.con.consess{c}.fcon.weights = [eye(phDSGN.nr_timebins) zeros(phDSGN.nr_timebins,size(noise_regs{1},2)) -eye(phDSGN.nr_timebins)];
        matlabbatch{3}.spm.stats.con.consess{c}.fcon.sessrep = 'none';
        c = c+1;
        matlabbatch{3}.spm.stats.con.consess{c}.fcon.name = [sess_subj{1} '_' sess_subj{2} '_interaction_effect'];
        matlabbatch{3}.spm.stats.con.consess{c}.fcon.weights = [diff(eye(phDSGN.nr_timebins)) zeros((phDSGN.nr_timebins-1),size(noise_regs{1},2)) -diff(eye(phDSGN.nr_timebins))];
        matlabbatch{3}.spm.stats.con.consess{c}.fcon.sessrep = 'none';
        
        for con = (b+1):c
            phDSGN.contrastnames{con} = matlabbatch{3}.spm.stats.con.consess{con}.fcon.name;
            phDSGN.contrastweights{con} = matlabbatch{3}.spm.stats.con.consess{con}.fcon.weights;
        end
    
    elseif n_contrasts == 3
        
        matlabbatch{3}.spm.stats.con.consess{c}.fcon.name = [sess_subj{1} '_main_effect'];
        matlabbatch{3}.spm.stats.con.consess{c}.fcon.weights = eye(phDSGN.nr_timebins);
        matlabbatch{3}.spm.stats.con.consess{c}.fcon.sessrep = 'none';
        c = c+1;
        matlabbatch{3}.spm.stats.con.consess{c}.fcon.name = [sess_subj{1} '_interaction_effect'];
        matlabbatch{3}.spm.stats.con.consess{c}.fcon.weights = diff(eye(phDSGN.nr_timebins));
        matlabbatch{3}.spm.stats.con.consess{c}.fcon.sessrep = 'none';
        c = c+1;
        matlabbatch{3}.spm.stats.con.consess{c}.fcon.name = [sess_subj{2} '_main_effect'];
        matlabbatch{3}.spm.stats.con.consess{c}.fcon.weights = [zeros(phDSGN.nr_timebins) zeros(phDSGN.nr_timebins,size(noise_regs{1},2)) eye(phDSGN.nr_timebins)];
        matlabbatch{3}.spm.stats.con.consess{c}.fcon.sessrep = 'none';
        c = c+1;
        matlabbatch{3}.spm.stats.con.consess{c}.fcon.name = [sess_subj{2} '_interaction_effect'];
        matlabbatch{3}.spm.stats.con.consess{c}.fcon.weights = [zeros(phDSGN.nr_timebins-1,phDSGN.nr_timebins) zeros(phDSGN.nr_timebins-1,size(noise_regs{1},2)) diff(eye(phDSGN.nr_timebins))];
        matlabbatch{3}.spm.stats.con.consess{c}.fcon.sessrep = 'none';
        c = c+1;
        matlabbatch{3}.spm.stats.con.consess{c}.fcon.name = [sess_subj{3} '_main_effect'];
        matlabbatch{3}.spm.stats.con.consess{c}.fcon.weights = [zeros(phDSGN.nr_timebins) zeros(phDSGN.nr_timebins,size(noise_regs{1},2)) zeros(phDSGN.nr_timebins) zeros(phDSGN.nr_timebins,size(noise_regs{2},2)) eye(phDSGN.nr_timebins)];
        matlabbatch{3}.spm.stats.con.consess{c}.fcon.sessrep = 'none';
        c = c+1;
        matlabbatch{3}.spm.stats.con.consess{c}.fcon.name = [sess_subj{3} '_interaction_effect'];
        matlabbatch{3}.spm.stats.con.consess{c}.fcon.weights = [zeros(phDSGN.nr_timebins-1,phDSGN.nr_timebins) zeros(phDSGN.nr_timebins-1,size(noise_regs{1},2)) zeros(phDSGN.nr_timebins-1,phDSGN.nr_timebins) zeros(phDSGN.nr_timebins-1,size(noise_regs{2},2)) diff(eye(phDSGN.nr_timebins))];
        matlabbatch{3}.spm.stats.con.consess{c}.fcon.sessrep = 'none';
        c = c+1;
        matlabbatch{3}.spm.stats.con.consess{c}.fcon.name = [sess_subj{1} '_' sess_subj{2} '_main_effect'];
        matlabbatch{3}.spm.stats.con.consess{c}.fcon.weights = [eye(phDSGN.nr_timebins) zeros(phDSGN.nr_timebins,size(noise_regs{1},2)) -eye(phDSGN.nr_timebins)];
        matlabbatch{3}.spm.stats.con.consess{c}.fcon.sessrep = 'none';
        c = c+1;
        matlabbatch{3}.spm.stats.con.consess{c}.fcon.name = [sess_subj{1} '_' sess_subj{2} '_interaction_effect'];
        matlabbatch{3}.spm.stats.con.consess{c}.fcon.weights = [diff(eye(phDSGN.nr_timebins)) zeros((phDSGN.nr_timebins-1),size(noise_regs{1},2)) -diff(eye(phDSGN.nr_timebins))];
        matlabbatch{3}.spm.stats.con.consess{c}.fcon.sessrep = 'none';
        c = c+1;
        matlabbatch{3}.spm.stats.con.consess{c}.fcon.name = [sess_subj{1} '_' sess_subj{3} '_main_effect'];
        matlabbatch{3}.spm.stats.con.consess{c}.fcon.weights = [eye(phDSGN.nr_timebins) zeros(phDSGN.nr_timebins,size(noise_regs{1},2)) zeros(phDSGN.nr_timebins) zeros(phDSGN.nr_timebins,size(noise_regs{2},2)) -eye(phDSGN.nr_timebins)];
        matlabbatch{3}.spm.stats.con.consess{c}.fcon.sessrep = 'none';
        c = c+1;
        matlabbatch{3}.spm.stats.con.consess{c}.fcon.name = [sess_subj{1} '_' sess_subj{3} '_interaction_effect'];
        matlabbatch{3}.spm.stats.con.consess{c}.fcon.weights = [diff(eye(phDSGN.nr_timebins)) zeros(phDSGN.nr_timebins-1,size(noise_regs{1},2)) zeros(phDSGN.nr_timebins-1,phDSGN.nr_timebins) zeros(phDSGN.nr_timebins-1,size(noise_regs{2},2)) -diff(eye(phDSGN.nr_timebins))];
        matlabbatch{3}.spm.stats.con.consess{c}.fcon.sessrep = 'none';
        c = c+1;
        matlabbatch{3}.spm.stats.con.consess{c}.fcon.name = [sess_subj{2} '_' sess_subj{3} '_main_effect'];
        matlabbatch{3}.spm.stats.con.consess{c}.fcon.weights = [zeros(phDSGN.nr_timebins) zeros(phDSGN.nr_timebins,size(noise_regs{1},2)) eye(phDSGN.nr_timebins) zeros(phDSGN.nr_timebins,size(noise_regs{2},2)) -eye(phDSGN.nr_timebins)];
        matlabbatch{3}.spm.stats.con.consess{c}.fcon.sessrep = 'none';
        c = c+1;
        matlabbatch{3}.spm.stats.con.consess{c}.fcon.name = [sess_subj{2} '_' sess_subj{3} '_interaction_effect'];
        matlabbatch{3}.spm.stats.con.consess{c}.fcon.weights = [zeros(phDSGN.nr_timebins-1,phDSGN.nr_timebins) zeros(phDSGN.nr_timebins-1,size(noise_regs{1},2)) diff(eye(phDSGN.nr_timebins)) zeros(phDSGN.nr_timebins-1,size(noise_regs{2},2)) -diff(eye(phDSGN.nr_timebins))];
        matlabbatch{3}.spm.stats.con.consess{c}.fcon.sessrep = 'none';
        
        for con = (b+1):c
            phDSGN.contrastnames{con} = matlabbatch{3}.spm.stats.con.consess{con}.fcon.name;
            phDSGN.contrastweights{con} = matlabbatch{3}.spm.stats.con.consess{con}.fcon.weights;
        end
        
    elseif n_contrasts > 3
        
        error('\n More than 3 sessions not yet implemented in script, please customize before proceeding\n');
        
    end
    
    matlabbatch{3}.spm.stats.con.delete = 0;
   
   
    % SAVE BATCHES AND RUN
    %---------------------
    
    cd(subfirstlevelmodeldir);
    
    eval(['save design_first_level_' derivsubjs{sub} '.mat matlabbatch']);

    spm_jobman('initcfg');
    spm_jobman('run',matlabbatch);
    
    cd(rootdir);
    
    
    % ADD CON IMAGES NAMES TO phDSGN STRUCTURE AND SAVE
    %--------------------------------------------------
    
    tconimgs = dir(fullfile(subfirstlevelmodeldir,'con*.nii'));
    fconimgs = dir(fullfile(subfirstlevelmodeldir,'ess*.nii'));
    conimgs = [tconimgs;fconimgs];
    
    phDSGN.contrasts = cell(1,size(conimgs,1));
    
    for conimg = 1:size(conimgs,1)
        phDSGN.contrasts{conimg} = fullfile(conimgs(conimg).folder,conimgs(conimg).name);
    end
    
    savefilename = fullfile(subfirstlevelmodeldir, 'DSGN.mat');
    save(savefilename, 'phDSGN', '-v7.3');

    
end