%% LaBGAScore_first_s1_prep_firstlvl.m
%
% SCRIPT NEEDS TO BE RUN FROM THE ROOT DIRECTORY OF YOUR DATASET
% 
% This script is an improved and integrated version of
% https://github.com/labgas/proj-whiplash/blob/main/prep/WAD_prep_s3_extract_confound_reg_fMRIprep.m
% https://github.com/labgas/proj-whiplash/blob/main/firstlevel/model_1_pictures/WAD_first_m1_s1_spm_prep_firstlvl_models.m
%
% It will extract noise regressors from fMRIprep output, including
% a) global CSF signal
% b) 24 head motion parameters (six directions, derivatives, and squared
% values)
% c) dummy spike regressors
% and save them as noise_regs.txt files, one per run
%
% It will also read events.tsv files and save them as onsets.mat files, one
% per run
%
% These will be saved in a correct folder structure for first level
% analysis with CANlab tools
% 
% DEPENDENCIES
% a) CANlabCore Github repo on your Matlab path
% b) SPM12 on your Matlab path
% c) BIDS data organization
%
% INPUTS 
% confound_regressor.tsv files from fMRIprep output
% raw func images (if spike_def = 'CANlab' - see below)
%
% OUTPUT
% noise_regs & onsets files that can be loaded into 
% - CANlab DSGN structure (RECOMMENDED) using LaBGAS_get_firstlvl_dsgn_obj.m
% - directly into SPM first level batch (https://github.com/labgas/proj-sert-fmri/blob/main/LaBGAS_first_level_batch_fMRIprep_conf.m)
%
% SPIKE_DEF (NOT CASE SENSITIVE)
% 'fMRIprep' use spike regressors based on a combination of DVARS and FD thresholds 
% set in fMRIprep arguments --fd-spike-threshold and --dvars-spike-threshold
%
% 'CANlab' use spike regressors based on CANlab's spike detection algorithm
% (Mahalanobis distance)(function scn_session_spike_id) and DVARS
% cfr make_nuisance_covs_from_fmriprep_output.m script in CANlab's
% CanLabScripts Github repo 
% https://github.com/canlab/CanlabScripts/tree/master/Scripts/Preprocessing
%
% DVARS_THRESHOLD
% set the threshold for standardized dvars to define a spike
% only used if spike_def = CANlab, otherwise set in fmriprep
% --dvars-spike-threshold command
% CANlab default is 3, but this is rather lenient
%
% OMIT_SPIKE_TRIALS
% 'no' do not remove onsets of pain trials coinciding with a spike
% 'yes' do remove - THIS IS NOT RECOMMENDED, WE PREFER TO DO THIS LATER
% BASED ON VIFS IN SINGLE TRIAL FIRST LEVEL ANALYSIS
%
% SPIKE_ADDITIONAL_VOLS
% set how many volumes after the spike you want to additionally regress out
% be careful for task-based data since this quite aggressive approach is
% mostly based on rs-fMRI, and beware of omitting too many volumes as well
% as creating missingness not at random - THIS IS NOT RECOMMENDED
%
% SPIKES_PERCENT_THRESHOLD
% set the maximum number of spikes (% of total volumes expressed as 0-1) you want to
% tolerate
%
% NOTE
% This script is NOT extensively tested for complex missing run patterns,
% unlike the scripts for proj-emosymp, as Iris' data hardly have any
% missing runs!
%__________________________________________________________________________
%
% authors: Lukas Van Oudenhove
% date:   March, 2022
%
%__________________________________________________________________________
% @(#)% LaBGAScore_first_s1_prep_firstlvl.m         v1.0        
% last modified: 2022/03/04


%% SET OPTIONS
%--------------------------------------------------------------------------

% MANDATORY
spike_def = 'fMRIprep';
omit_spike_trials = 'no';
spikes_percent_threshold=0.15;

% ONLY NEEDED IF SPIKE_DEF = CANlab
dvars_threshold = 2;
spike_additional_vols=0;

% THIS CHOICE OF OPTIONS CAN BE CONSIDERED LABGAS DEFAULTS, BUT MAY BE
% STUDY SPECIFIC, SO DISCUSS WITH LUKAS IF IN DOUBT!


%% DEFINE AND CREATE DIRECTORIES AND RUNS
%--------------------------------------------------------------------------

% load standard BIDS directory structure from root dir

LaBGAScore_prep_s0_define_directories;

% create first level directory

firstleveldir = fullfile(rootdir,'firstlevel');
if ~exist(firstleveldir,'dir')
    mkdir(firstleveldir);
end

% define runs

runs = {'run-1';'run-2';'run-3';'run-4';'run-5';'run-6'};


%% CREATE CANLAB DSGN STRUCTURE
    
% INPUT
    % required fields
    DSGN.metadata = "proj-erythritol_4a first level analysis model 1, i.e. modeling 4 conditions for high, moderate, neutral fear as long events (= duration of solution in mouth)"; % field for annotation with study info, or whatever you like
    DSGN.modeldir = '/data/test_scripts/firstlevel/model_1_basic'; % directory where you want to write first level results for this model
        if ~exist(DSGN.modeldir,'dir')
            mkdir(DSGN.modeldir);
        end
    DSGN.subjects = derivsubjdirs'; % sets up empty cell array field for subjects in structure array DSGN
    DSGN.funcnames = {'/func/run-1/s6*.nii',...
        '/func/run-2/s6*.nii',...
        '/func/run-3/s6*.nii',...
        '/func/run-4/s6*.nii',...
        '/func/run-5/s6*.nii',...
        '/func/run-6/s6*.nii'}; % cell array (one cell per session) of paths to functional files, relative to absolute path specific in DSGN.subjects
   
    % optional fields
    DSGN.concatenation = {}; % default: none; cell array of arrays of runs to concatenate; see documentation for when to concatenate, and how it works exactly
    DSGN.allowmissingfunc = true; % default: false; true will prevent erroring out when functional file is missing for at least one run is missing for at least one subject
    DSGN.customrunintercepts = {}; % default: none; will only work if DSGN.concatenation is specified; cell array of vectors specifying custom intercepts
    
% PARAMETERS
    DSGN.tr = 1.8; % repetition time (TR) in seconds
    DSGN.hpf = 128; % high pass filter in seconds; SPM default is 128, CANlab default is 180 since the brain response to pain stimuli last long and variance may be lost at shorter lengths, use scn_spm_design_check output, and review the SPM.mat in spm for diagnostics; 
    % STUDY-SPECIFIC: in this study, with rather short non-pain events and
    % relatively short ITI, we stick with the SPM default
    DSGN.fmri_t = 30; % microtime resolution - t=number of slices since we did slice timing; spm (and CANlab) default 16 can be kept for multiband w/o slice timing; TO BE CHECKED SINCE WE HAVE MULTIBAND WITH SLICE TIMING
    DSGN.fmri_t0 = 15; % microtime onset - reference slice used in slice timing correction; spm (and CANlab) default 1 can be kept for multiband w/o slice timing
    
% MODELING
    % required fields
    
    % cell array (one cell per session (i.e. run in our case)) of cell
    % arrays (one cell per condition) of MAT-file names, in fixed order:
    % all conditions of interest first, conditions of no interest last
    % if only one session is specified, it will be applied to all sessions
    % TO BE CHECKED WHETHER THIS WORKS
    c=0;
    c=c+1;DSGN.conditions{c}={'sucrose' 'erythritol' 'sucralose' 'water' 'rating' 'swallow_rinse'};
    
    % optional fields
%     DSGN.pmods = {{}}; % cell array (one cell per session) of cell arrays (one cell per condition) of cell arrays (one cell per modulator) of MAT-file names
%     DSGN.convolution; default hrf.derivs = [0 0]; structure specifying the convolution to use for conditions different fields required depending on convolution type; 
%     DSGN.ar1 = false; % autoregressive AR(1) to model serial correlations; SPM default is true, CANlab default is false, Tor recommends turning autocorrelation off, because this algorithm pools across the whole brain, and does not perform well in some situations; if you are performing a group analysis, the autocorrelation problem is not as concerning
    DSGN.notimemod = true; % default: false; if true, turn off time modulation of conditions, i.e. when you do not expect linear trends over time
%     DSGN.singletrials = {{}}; % a cell array (1 cell per session) of cell arrays (1 cell per condition) of (corresponding to DSGN.conditions) of true/false values indicating whether to convert specified condition to set of single trial conditions
%     DSGN.singletrialsall = false; % default: false; if true, set DSGN.singletrials to true for all conditions
    DSGN.modelingfilesdir = 'model_1_basic'; % name of subfolder which will be created within directory containing functional files where .mat files containing fields of DSGN structure will be saved
%     DSGN.allowemptycond = false; % default:false; if true, allow empty conditions
%     DSGN.allowmissingcondfiles = false; % default:false; if true, throw warning instead of error when no file(s) are found corresponding to a MAT-file name/wildcard
    DSGN.multireg = 'noise_regs'; % specify name for matfile with noise parameters you want to save
    
    % CONTRASTS
    % required fields
    % cell array (one cell per contrast) of contrast definitions
    c=0;
    c=c+1;DSGN.contrasts{c} = {{'sucrose'}}; % CON_0001
    DSGN.contrastnames{c} = 'sucrose'; % not needed strictly, because this will be automatically generated for standard contrasts like this (sucrose versus baseline)
    DSGN.contrastweights{c} = [1];
    c=c+1;DSGN.contrasts{c} = {{'erythritol'}}; % CON_0002
    DSGN.contrastnames{c} = 'erythritol'; % not needed strictly, because this will be automatically generated for standard contrasts like this (erythritol versus baseline)
    DSGN.contrastweights{c} = [1];
    c=c+1;DSGN.contrasts{c} = {{'sucralose'}}; % CON_0003
    DSGN.contrastnames{c} = 'sucralose'; % not needed strictly, because this will be automatically generated for standard contrasts like this (sucralose versus baseline)
    DSGN.contrastweights{c} = [1];
    c=c+1;DSGN.contrasts{c} = {{'water'}}; % CON_0004
    DSGN.contrastnames{c} = 'water'; % not needed strictly, because this will be automatically generated for standard contrasts like this (sucralose versus baseline)
    DSGN.contrastweights{c} = [1];
    c=c+1;DSGN.contrasts{c} = {{'sucrose'} {'water'}}; % CON_0005
    DSGN.contrastnames{c} = 'sucrose_vs_water'; % not needed strictly, because this will be automatically generated for standard contrasts like this (sucrose versus water)
    DSGN.contrastweights{c} = [1 -1];
    c=c+1;DSGN.contrasts{c} = {{'erythritol'} {'water'}}; % CON_0006
    DSGN.contrastnames{c} = 'erythritol_vs_water'; % not needed strictly, because this will be automatically generated for standard contrasts like this (erythritol versus water)
    DSGN.contrastweights{c} = [1 -1];
    c=c+1;DSGN.contrasts{c} = {{'sucralose'} {'water'}}; % CON_0007
    DSGN.contrastnames{c} = 'sucralose_vs_water'; % not needed strictly, because this will be automatically generated for standard contrasts like this (sucralose versus water)
    DSGN.contrastweights{c} = [1 -1];
    c=c+1;DSGN.contrasts{c} = {{'sucrose'} {'sucralose'}}; % CON_0008
    DSGN.contrastnames{c} = 'sucrose_vs_sucralose'; % not needed strictly, because this will be automatically generated for standard contrasts like this (sucrose versus sucralose)
    DSGN.contrastweights{c} = [1 -1];
    c=c+1;DSGN.contrasts{c} = {{'sucrose'} {'erythritol'}}; % CON_0009
    DSGN.contrastnames{c} = 'sucrose_vs_erythritol'; % not needed strictly, because this will be automatically generated for standard contrasts like this (sucrose versus erythritol)
    DSGN.contrastweights{c} = [1 -1];
    c=c+1;DSGN.contrasts{c} = {{'erythritol'} {'sucralose'}}; % CON_0010
    DSGN.contrastnames{c} = 'erythritol_vs_sucralose'; % not needed strictly, because this will be automatically generated for standard contrasts like this (erythritol versus sucralose)
    DSGN.contrastweights{c} = [1 -1];

%--------------------------------------------------------------------------
% END OF STUDY-SPECIFIC CODE, BELOW SHOULD WORK OUT OF THE BOX
%--------------------------------------------------------------------------

%% CREATE MODELDIR IN FIRSTLEVELDIR

% create firstmodeldir
firstmodeldir = fullfile(firstleveldir,DSGN.modelingfilesdir);
if ~exist(firstmodeldir,'dir')
    mkdir(firstmodeldir);
end

cd (firstmodeldir);

% write subjectdirs in firstlevelmodeldir
sm=@(x)mkdir(x); % defines mkdir as an anonymous function sm
cellfun(sm,derivsubjs);


%% LOOP OVER SUBJECTS
%--------------------------------------------------------------------------

for sub=1:size(derivsubjs,1)
    
    %% WRITE RUNDIRS & DEFINE SUBJECT LEVEL FILENAMES
    
    subjderivdir = fullfile(derivsubjdirs{sub},'func');
    cd(subjderivdir);
    cellfun(sm,runs);
    
    
    rawimgs = dir(fullfile(BIDSsubjdirs{sub},'func/*bold.nii.gz'));
    rawimgs = {rawimgs(:).name}';
    if strcmpi(spike_def,'CANlab')==1 % we only need raw images in case spike_def = CANlab
        for rawimg = 1:size(rawimgs,1) % raw images are needed when spike_def = CANlab, which calls a function that is incompatible with .nii.gz, hence we unzip
            gunzip(rawimgs{rawimg});
        end
        rawimgs = dir(fullfile(BIDSsubjdirs{sub},'func/*bold.nii'));
        rawimgs = {rawimgs(:).name}';
        rawidx = ~contains(rawimgs,'rest'); % omit resting state scan if it exists
        rawimgs = {rawimgs{rawidx}}';
    else
        continue
    end
    
    preprocimgs = dir(fullfile(subjderivdir,'s6-*.nii.gz'));
    preprocimgs = {preprocimgs(:).name}';
        for preprocimg = 1:size(preprocimgs,1) % raw images are needed when spike_def = CANlab, which calls a function that is incompatible with .nii.gz, hence we unzip
                gunzip(preprocimgs{preprocimg});
                delete(preprocimgs{preprocimg});
        end
    preprocimgs = dir(fullfile(subjderivdir,'s6-*.nii'));
    preprocimgs = {preprocimgs(:).name}';
    preprocidx = ~contains(preprocimgs,'rest'); % omit resting state scan if it exists
    preprocimgs = {preprocimgs{preprocidx}}';
    
    fmriprep_noisefiles = dir(fullfile(subjderivdir,'*desc-confounds_timeseries.tsv'));
    fmriprep_noisefiles = {fmriprep_noisefiles(:).name}';
    noiseidx = ~contains(fmriprep_noisefiles,'rest'); % omit resting state scan if it exists
    fmriprep_noisefiles = {fmriprep_noisefiles{noiseidx}}';
    
    runnames = char(fmriprep_noisefiles);
    runnames = runnames(:,1:29); % this is study-specific, depends on length of subjectname and taskname - could be adapted based on regexp, but easy enough to adapt per study
        if ~isequal(size(rawimgs,1),size(preprocimgs,1),size(fmriprep_noisefiles,1)) 
            error('numbers of raw images, preprocessed images, and noise files do not match for %s, please check before proceeding',derivsubjs{sub});
        else
            disp('numbers of raw images, preprocessed images, and noise files match for %s, good to go',derivsubjs{sub});
        end

    %% CALCULATE AND/OR EXTRACT CONFOUND REGRESSORS AND EVENTS FILES, AND WRITE TO THE CORRECT FILE/FOLDER STRUCTURE
    
    for run=1:size(fmriprep_noisefiles,1)
        
        % define subdir for this run
        rundir = fullfile(subjderivdir,runs{run});
        
        if contains(fmriprep_noisefiles{run},runs{run}) % check whether data for run are not missing/excluded
        
            % move preprocessed image and fmriprep noisefile into rundir
            movefile(fullfile(subjderivdir,fmriprep_noisefiles{run}),fullfile(rundir,fmriprep_noisefiles{run}));
            movefile(fullfile(subjderivdir,preprocimgs{run}),fullfile(rundir,preprocimgs{run}));

            % CONFOUND REGRESSOR FILES
            % load confound regressor file generated by fMRIprep into Matlab table
            % variable
            Rfull = readtable(fullfile(rundir,fmriprep_noisefiles{run}),'TreatAsEmpty','n/a','FileType', 'text', 'Delimiter', 'tab');

            % replace NaNs in first row with Os
            wh_replace = ismissing(Rfull(1,:));
                if any(wh_replace)
                    Rfull{1, wh_replace} = zeros(1, sum(wh_replace)); % make array of zeros of the right size
                end

            % calculate and extract confound regressors
                if strcmpi(spike_def,'fMRIprep')==1 % switch would probably make more sense in this case, but this works too!

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
                    % 'spike_additional_vols' remains unspecified, everything will proceed as
                    % it did before, meaning spikes will be identified and flagged in the
                    % creation of nuisance regressors without considering the following TRs
                    % Add them if user requested, for both nuisance_covs and dvars_spikes_regs
                        if exist('spike_additional_vols')
                            additional_spikes_regs = zeros(height(Rfull),size(spikes,2)*spike_additional_vols);
                                % This loop will create a separate column with ones in each row (TR) 
                                % we would like to consider a nuisance regressor
                                for i = 1:size(spikes,2) 
                                    additional_spikes_regs(spikes(i)+1 : spikes(i)+spike_additional_vols,(i*spike_additional_vols-(spike_additional_vols-1)):(i*spike_additional_vols)) = eye(spike_additional_vols);
                                end
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
                        for i = 1:length(duplicate_rows) % This loop sets duplicate values to zero; drops them later (to keep indices the same during the loop)
                            [~,curr_cols] = find(Rfull{duplicate_rows(i),:}==1);
                            Rfull{duplicate_rows(i), curr_cols(2:end)} = 0;
                        end
                    Rfull = Rfull(1:height(Rfull), any(table2array(Rfull)));

                elseif strcmpi(spike_def,'CANlab')==1

                    % define raw image file
                    raw_img_fname = rawimgs{run};

                    % add in canlab spike detection (Mahalanobis distance)
                    [g, mahal_spikes, gtrim, mahal_spikes_regs, snr] = scn_session_spike_id(fullfile(subjrawdir,raw_img_fname), 'doplot', 0); % CANlab function needs to be on your Matlab path
                    delete('*.img'); % delete implicit mask .hdr/.img files generated by the CANlab function on the line above, since we don't need/use them
                    delete('*.hdr');
                    mahal_spikes_regs(:,1) = []; %drop gtrim which is the global signal
                    Rfull(:,contains(Rfull.Properties.VariableNames,'motion_outlier'))=[]; % drop fmriprep motion outliers since we do not use them when spike_def = CANlab, and they cause redundancies
                    Rfull = [Rfull array2table(mahal_spikes_regs)];

                    % add in dvars spike regressors that are non-redundant with mahal spikes
                    dvarsZ = [0; zscore(Rfull.dvars(2:end))]; % first element of dvars always = 0, drop it from zscoring and set it to Z=0
                    dvars_spikes = find(dvarsZ > dvars_threshold);
                    same = ismember(dvars_spikes,mahal_spikes);
                    dvars_spikes(same) = []; % drop the redundant ones
                    dvars_spikes_regs = zeros(height(Rfull),size(dvars_spikes,1));
                        for i=1:size(dvars_spikes,1)
                            dvars_spikes_regs(dvars_spikes(i),i) = 1;
                        end
                    Rfull = [Rfull array2table(dvars_spikes_regs)];

                    % flag user-specified number of volumes after each spike
                    % Motion can create artifacts lasting longer than the single image we
                    % usually account for using spike id scripts. we're also going to flag the
                    % following TRs, the number of which is defined by the user. If
                    % 'spike_additional_vols' remains unspecified, everything will proceed as
                    % it did before, meaning spikes will be identified and flagged in the
                    % creation of nuisance regressors without considering the following TRs
                    % Add them if user requested, for both nuisance_covs and dvars_spikes_regs
                        if exist('spike_additional_vols')
                            % concatenate generated spike and DVARS regs. We
                            % would like to flag subsequent TR's with respect to both of these
                            % measures.
                            spikes = [mahal_spikes;dvars_spikes];
                            additional_spikes_regs = zeros(size(mahal_spikes_regs,1),size(spikes,1)*spike_additional_vols);
                                % This loop will create a separate column with ones in each row (TR) 
                                % we would like to consider a nuisance regressor
                                % Performs this function for spikes and DVARS. 
                                for i = 1:size(spikes,1) 
                                    additional_spikes_regs(spikes(i)+1 : spikes(i)+spike_additional_vols,(i*spike_additional_vols-(spike_additional_vols-1)):(i*spike_additional_vols)) = eye(spike_additional_vols);
                                end
                            % if any spikes went beyond the end, trim it down
                            additional_spikes_regs = additional_spikes_regs(1:height(Rfull),:);
                            % add the additional spikes to the larger matrix
                            Rfull = [Rfull array2table(additional_spikes_regs)];
                        end

                    % remove redundant spike regressors
                    regs = Rfull.Properties.VariableNames;
                    spike_cols = contains(regs,'mahal_spikes'); 
                    dvars_cols = contains(regs,'dvars_spikes'); 
                    additional_spike_cols = contains(regs,'additional_spikes'); 

                    [duplicate_rows, ~] = find(sum(Rfull{:, spike_cols | dvars_cols | additional_spike_cols}, 2)>1);
                        for i = 1:size(duplicate_rows,1) %This loop sets duplicate values to zero; drops them later (to keep indices the same during the loop)
                            [~,curr_cols] = find(Rfull{duplicate_rows(i),:}==1);
                            Rfull{duplicate_rows(i), curr_cols(2:end)} = 0;
                        end
                    Rfull = Rfull(1:size(mahal_spikes_regs,1), any(table2array(Rfull)));
                else
                    error('invalid spike_def option')
                end

            % Select confound and spike regressors to return for use in GLM 
            regs = Rfull.Properties.VariableNames;
            motion_cols = contains(regs,'rot') | contains(regs,'trans');
            spike_cols = contains(regs,'mahal_spikes') | contains(regs,'motion_outlier'); 
            dvars_cols = contains(regs,'dvars_spikes'); 
            additional_spike_cols = contains(regs,'additional_spikes'); 
            R = Rfull(:,motion_cols | spike_cols | dvars_cols | additional_spike_cols);
            R.csf = Rfull.csf;
            Rspikes=Rfull(:,spike_cols | dvars_cols | additional_spike_cols);
            Rspikes.spikes=sum(Rspikes{:,1:end},2);
            volume_idx = [1:height(Rfull)]; 
            spikes = volume_idx(Rspikes.spikes==1)';

            % compute and output how many spikes total
            n_spike_regs = sum(dvars_cols | spike_cols | additional_spike_cols);
            n_spike_regs_percent = n_spike_regs / height(Rfull);

            % print warning if #volumes identified as spikes exceeds
            % user-defined threshold
                if n_spike_regs_percent > spikes_percent_threshold
                    warning('number of volumes identified as spikes exceeds threshold in %s',runnames(run,:))
                end

            % save confound regressors as matrix named R for use in
            % SPM/CANlab GLM model tools
            R=table2array(R);
            filename_noise_regs = fullfile(rundir,DSGN.multireg);
            save(filename_noise_regs,'R');
            
            clear R* filename_events

            % EVENTS FILES
            % read events.tsv files with onsets, durations, and trial type
            eventsfiles = dir(fullfile(subjrawdir,strcat(runnames(run,:),'*events.tsv')));
            eventsfiles = {eventsfiles(:).name}';
            % loop over models
            for model = 1:size(models,1)
                O=readtable(fullfile(subjrawdir,eventsfiles{model}),'FileType', 'text', 'Delimiter', 'tab');
                O.trial_type = categorical(O.trial_type);
                % omit trials that coincide with spikes if that option is chosen
                    if strcmpi(omit_spike_trials,'yes')==1
                        same=ismember(O.onset,spikes); % identify trials for which onset coincides with spike
                        O(same,:)=[]; % get rid of trials coinciding with spikes
                    elseif strcmpi(omit_spike_trials,'no')==1
                    else
                        error('invalid omit_spike_trials option')
                    end
                % save events file as .mat file
                filename_events = strcat(rundir,'\onsets_',runnames(run,:),'_',models{model});
                save(filename_events,'O');
                clear O filename_events
            end % for loop models

            clear onsetfiles
            
        else
            continue
            
        end % if loop checking if data for run is not missing
        
    end % for loop runs
    
end % for loop subjects 