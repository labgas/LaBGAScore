%% LaBGAScore_firstlevel_s1b_fit_phMRI_model.m
%
%
% *USAGE*
%
% This scripts performs first-level model specification, estimation, and
% contrast definition and estimation for phMRI designs with any number of
% conditions (one session per condition).
%
% The noise model per session is the six realignment parameters and their
% first derivatives, z-scored (12 columns), optionally their squares
% (24 columns), the fMRIprep motion_outlier spike regressors if
% phDSGN.model_spikes is true, and the global CSF signal.
%
% *OPTIONS*
%
% For LaBGAS_options fields (subjs2analyze, mandatory.spikes_percent_threshold,
% movement_reg_quadratic, spikes.spike_additional_vols), see
% LaBGAScore_firstlevel_s1_options_dsgn_struct.m
%
% phDSGN fields (this script's own struct, all STUDY-SPECIFIC unless noted):
%
% * modelingfilesdir    name of subfolder within firstlevel dir where model output is saved
% * tr                  repetition time (TR) in seconds
% * t                   microtime resolution, set to nr of slices if slice timing was performed
% * t0                  microtime onset - reference slice used in slice timing correction
% * hpf                 high pass filter length in seconds; very long in case of phMRI
% * multireg            name for the noise-regressor matfile written per session
% * nr_dyns             nr of dynamics (volumes) per run
% * start_dyn           dynamic (volume) number where administration starts
% * timebin_length      length of each post-administration timebin, in seconds
% * sessions            cell array of condition names, in order of corresponding sessions
% * model_spikes        true/false, add fMRIprep motion_outlier columns as spike regressors
%
% *DEPENDENCIES*
%
% SPM12 on Matlab path (spm_jobman, cfg_dep, spm_select)
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
% LaBGAScore_firstlevel_s1b_fit_phMRI_model.m         v1.3
%
% last modified: 2026/10/30
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

% noise model
phDSGN.model_spikes = true;    % add fMRIprep motion_outlier columns as spike regressors; false gives 24hmp (or 12hmp) + csf only

% calculations, don't change
phDSGN.timebin_length_tr = phDSGN.timebin_length/phDSGN.tr;
phDSGN.nr_timebins = (phDSGN.nr_dyns-phDSGN.start_dyn)/phDSGN.timebin_length_tr; % nr of post-infusion timebins

    if mod(phDSGN.nr_timebins,1) ~= 0 || mod(phDSGN.timebin_length_tr,1) ~= 0
        error('\ntimebin_length must be a whole number of TRs and (nr_dyns - start_dyn) a whole number of timebins, please check phDSGN before proceeding\n');
    end


% LABGAS OPTIONS

LaBGAS_options = struct();

% subjects
LaBGAS_options.subjs2analyze = {}; % enter subjects separated by comma if you only want to analyze selected subjects e.g. {'sub-001','sub-002'}; leave empty to loop over all subjects

% required options
LaBGAS_options.mandatory.spikes_percent_threshold=0.15;
LaBGAS_options.movement_reg_quadratic = true; % change to false if you don't want to add quadratic terms for movement parameters and their first-order derivatives

% spike options, only used if phDSGN.model_spikes is true
LaBGAS_options.spikes.spike_additional_vols=0; % OPTIONAL, NOT RECOMMENDED TO TURN ON


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
        phDSGN.modeldir = firstlevelmodeldir; % make CANlab DSGN structure for the purpose of using it in second level analysis
    end
    
phDSGN.subjects = derivsubjdirs';

    if ~isempty(LaBGAS_options.subjs2analyze)
        missing_subjs = setdiff(LaBGAS_options.subjs2analyze,derivsubjs);
            if ~isempty(missing_subjs)
                error('\nsubject(s) %s in LaBGAS_options.subjs2analyze not found in %s, please check before proceeding\n',strjoin(missing_subjs,', '),derivdir);
            end
        derivsubjs = LaBGAS_options.subjs2analyze(:);
    end

    
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
    phDSGN.conditions = {};
    phDSGN.contrastnames = {}; % reset per subject, so nothing is left over from a previous subject with more sessions
    phDSGN.contrastweights = {};
    phDSGN.contrasts = {};

        for sess = 1:size(phDSGN.sessions,2)

            derivsubjsesdirs{sess} = fullfile(derivdir,derivsubjs{sub},['ses-' num2str(sess)],'func'); % built from the subject NAME, so a LaBGAS_options.subjs2analyze subset stays aligned
            derivsubjsesmodeldirs{sess} = fullfile(derivsubjsesdirs{sess},phDSGN.modelingfilesdir);
                if exist(derivsubjsesdirs{sess},'dir')
                    if ~exist(derivsubjsesmodeldirs{sess},'dir')
                        mkdir(derivsubjsesmodeldirs{sess});
                    end
                end
                
            derivimg = dir(fullfile(derivsubjsesdirs{sess},'s6*.nii.gz'));
                if isempty(derivimg)
                    derivimg = dir(fullfile(derivsubjsesdirs{sess},'s6*.nii')); % smoothed image may only exist unzipped, e.g. after the .nii.gz was removed
                end
                if size(derivimg,1) > 1
                    error('\nMore than one image file found with filter s6*.nii(.gz) in dir %s, please check before proceeding\n', derivsubjsesdirs{sess});
                elseif size(derivimg,1) < 1
                    fprintf('\nWARNING: No image files found with filter s6*.nii(.gz) in dir %s, missing condition - proceeding\n', derivsubjsesdirs{sess});
                    continue
                else
                    derivimg_name = regexprep(derivimg.name,'\.gz$','');
                        if ~exist(fullfile(derivimg.folder,derivimg_name),'file')
                            gunzip(fullfile(derivimg.folder,derivimg.name));
                        end
                    derivimages{sess} = spm_select('ExtFPlist', derivimg.folder, ['^' regexptranslate('escape',derivimg_name) '$'], Inf);
                        if size(derivimages{sess},1) ~= phDSGN.nr_dyns
                            error('\n%d volumes found in %s, phDSGN.nr_dyns is %d; check that the image content is present (datalad get) and nr_dyns is correct before proceeding\n', size(derivimages{sess},1), fullfile(derivimg.folder,derivimg_name), phDSGN.nr_dyns);
                        end
                end

            fmriprep_noisefile = dir(fullfile(derivsubjsesdirs{sess},'*desc-confounds_timeseries.tsv'));
                if size(fmriprep_noisefile,1) > 1
                    error('\nMore than one noise file found with filter *desc-confounds_timeseries.tsv in dir %s, please check before proceeding\n', derivsubjsesdirs{sess});
                elseif size(fmriprep_noisefile,1) < 1
                    fprintf('\nWARNING: No noise files found with filter *desc-confounds_timeseries.tsv in dir %s, missing condition - proceeding\n', derivsubjsesdirs{sess});
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
                    if height(Rfull) ~= phDSGN.nr_dyns
                        error('\n%d rows found in %s, phDSGN.nr_dyns is %d, please check before proceeding\n', height(Rfull), fmriprep_noisefiles{sess}, phDSGN.nr_dyns);
                    end

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
                zscorefun = @(x) zscore(x); % not named zscore: a variable of that name shadows the function on every later call in this workspace
                Rmotionzscore = varfun(zscorefun,Rmotionderiv);
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
                    if phDSGN.model_spikes
                        R = [Rmotionfinal,Rspikes,Rcsf];

                        % compute and output how many spikes total
                        n_spike_regs = width(Rspikes);
                        n_spike_regs_percent = n_spike_regs / height(Rfull);
                        fprintf('\n%s ses-%d: %d motion, %d spike and 1 csf regressor(s)\n', derivsubjs{sub}, sess, width(Rmotionfinal), n_spike_regs);

                        % print warning if #volumes identified as spikes exceeds
                        % user-defined threshold
                            if n_spike_regs_percent > LaBGAS_options.mandatory.spikes_percent_threshold
                                warning('number of volumes identified as spikes (%1.2f) exceeds threshold %1.2f in %s ses-%d',n_spike_regs_percent,LaBGAS_options.mandatory.spikes_percent_threshold,derivsubjs{sub},sess)
                            end
                    else
                        R = [Rmotionfinal,Rcsf];
                        fprintf('\n%s ses-%d: %d motion and 1 csf regressor(s), spikes not modelled\n', derivsubjs{sub}, sess, width(Rmotionfinal));
                    end

                names = R.Properties.VariableNames;

                % save confound regressors as matrix named R for use in
                % SPM/CANlab GLM model tools
                R=table2array(R);
                noise_regs{sess} = R;

                % write confound regressors
                filename_noise_regs = fullfile(derivsubjsesmodeldirs{sess},phDSGN.multireg);
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
    noise_regs = noise_regs(~empty_derivimgs); % keep aligned with sess_subj, otherwise the contrast weights skip the wrong noise columns when an earlier session is missing
    
        for sess = 1:size(sess_subj,2)
            
            phDSGN.conditions{sess} = cell(1,phDSGN.nr_timebins);
            
            matlabbatch{1}.spm.stats.fmri_spec.sess(sess).scans = cellstr(derivimages{sess});
            
            for b = 1:phDSGN.nr_timebins
                
                matlabbatch{1}.spm.stats.fmri_spec.sess(sess).cond(b).name = [sess_subj{sess} '_bin_' num2str(b)];
                matlabbatch{1}.spm.stats.fmri_spec.sess(sess).cond(b).onset = phDSGN.timebin_length_tr * b + (phDSGN.start_dyn - phDSGN.timebin_length_tr);
                matlabbatch{1}.spm.stats.fmri_spec.sess(sess).cond(b).duration = phDSGN.timebin_length_tr;
                matlabbatch{1}.spm.stats.fmri_spec.sess(sess).cond(b).tmod = 0;
                matlabbatch{1}.spm.stats.fmri_spec.sess(sess).cond(b).pmod = struct('name', {}, 'param', {}, 'poly', {});
                matlabbatch{1}.spm.stats.fmri_spec.sess(sess).cond(b).orth = 1;
                
                phDSGN.conditions{sess}{b} = [sess_subj{sess} '_bin_' num2str(b)];
                
            end
            
            matlabbatch{1}.spm.stats.fmri_spec.sess(sess).hpf = phDSGN.hpf;
            matlabbatch{1}.spm.stats.fmri_spec.sess(sess).multi_reg = {fullfile(derivsubjsesmodeldirs{sess},[phDSGN.multireg '.mat'])};
            
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
    
    % Each session's columns in the design matrix are its nr_timebins
    % condition regressors followed by its own noise regressors, whose number
    % differs between sessions when spikes are modelled. sess_cols{s} holds
    % the condition columns of session s, so every weight is placed by
    % column index rather than by counting zeros.
    
    matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    
    nr_sess_subj = size(sess_subj,2);
    nr_cols_sess = phDSGN.nr_timebins + cellfun(@(x) size(x,2), noise_regs);
    first_col_sess = [0 cumsum(nr_cols_sess(1:end-1))];
    nr_cols = sum(nr_cols_sess); % session constants come after this and get weight 0
    sess_cols = cell(1,nr_sess_subj);
        for s = 1:nr_sess_subj
            sess_cols{s} = first_col_sess(s) + (1:phDSGN.nr_timebins);
        end
        
        if nr_sess_subj >= 2
            sess_pairs = nchoosek(1:nr_sess_subj,2); % (1,2), (1,3), (2,3) for three sessions, the order of the original hardcoded contrasts
        else
            sess_pairs = zeros(0,2);
            fprintf('\nWARNING: Only 1 session found for subject %s, writing within-session F-contrasts only\n', derivsubjs{sub});
        end
    
    con = 0;
    
    % T-CONTRASTS: session i vs session j, per timebin
    
    for p = 1:size(sess_pairs,1)
        
        i = sess_pairs(p,1);
        j = sess_pairs(p,2);
        
        for a = 1:phDSGN.nr_timebins
            
            con = con + 1;
            weights = zeros(1,nr_cols);
            weights(sess_cols{i}(a)) = 1;
            weights(sess_cols{j}(a)) = -1;
            
            matlabbatch{3}.spm.stats.con.consess{con}.tcon.name = [sess_subj{i} '_' sess_subj{j} '_bin_' num2str(a)];
            matlabbatch{3}.spm.stats.con.consess{con}.tcon.weights = weights;
            matlabbatch{3}.spm.stats.con.consess{con}.tcon.sessrep = 'none';
            
        end
        
    end
    
    % F-CONTRASTS: main effect and interaction (effect of timebin) per session,
    % then the same for each pair of sessions
    
    main_w = eye(phDSGN.nr_timebins);
    inter_w = diff(eye(phDSGN.nr_timebins));
    
    nr_fcons = 2*nr_sess_subj + 2*size(sess_pairs,1);
    fcon_names = cell(1,nr_fcons);
    fcon_weights = cell(1,nr_fcons);
    f = 0;
    
        for s = 1:nr_sess_subj
            
            weights = zeros(phDSGN.nr_timebins,nr_cols);
            weights(:,sess_cols{s}) = main_w;
            f = f + 1;
            fcon_names{f} = [sess_subj{s} '_main_effect'];
            fcon_weights{f} = weights;
            
            weights = zeros(phDSGN.nr_timebins-1,nr_cols);
            weights(:,sess_cols{s}) = inter_w;
            f = f + 1;
            fcon_names{f} = [sess_subj{s} '_interaction_effect'];
            fcon_weights{f} = weights;
            
        end
        
        for p = 1:size(sess_pairs,1)
            
            i = sess_pairs(p,1);
            j = sess_pairs(p,2);
            
            weights = zeros(phDSGN.nr_timebins,nr_cols);
            weights(:,sess_cols{i}) = main_w;
            weights(:,sess_cols{j}) = -main_w;
            f = f + 1;
            fcon_names{f} = [sess_subj{i} '_' sess_subj{j} '_main_effect'];
            fcon_weights{f} = weights;
            
            weights = zeros(phDSGN.nr_timebins-1,nr_cols);
            weights(:,sess_cols{i}) = inter_w;
            weights(:,sess_cols{j}) = -inter_w;
            f = f + 1;
            fcon_names{f} = [sess_subj{i} '_' sess_subj{j} '_interaction_effect'];
            fcon_weights{f} = weights;
            
        end
        
        for f = 1:numel(fcon_names)
            
            con = con + 1;
            
            matlabbatch{3}.spm.stats.con.consess{con}.fcon.name = fcon_names{f};
            matlabbatch{3}.spm.stats.con.consess{con}.fcon.weights = fcon_weights{f};
            matlabbatch{3}.spm.stats.con.consess{con}.fcon.sessrep = 'none';
            
        end
        
        for c = 1:con
            consess = matlabbatch{3}.spm.stats.con.consess{c};
            contype = fieldnames(consess);
            phDSGN.contrastnames{c} = consess.(contype{1}).name;
            phDSGN.contrastweights{c} = consess.(contype{1}).weights;
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
    
    % taken from SPM.xCon, so phDSGN.contrasts{c} is the image of
    % phDSGN.contrastnames{c} (con_*.nii for t-, ess_*.nii for F-contrasts)
    
    load(fullfile(subfirstlevelmodeldir,'SPM.mat'),'SPM');
    
    phDSGN.contrasts = cell(1,size(SPM.xCon,2));
    
    for conimg = 1:size(SPM.xCon,2)
        phDSGN.contrasts{conimg} = fullfile(subfirstlevelmodeldir,SPM.xCon(conimg).Vcon.fname);
    end
    
    clear SPM
    
    savefilename = fullfile(subfirstlevelmodeldir, 'DSGN.mat');
    save(savefilename, 'phDSGN', '-v7.3');

    
end
