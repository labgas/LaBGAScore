%% WAD_prep_s3_extract_confound_reg_fMRIprep
% 
% This script is an improved and integrated version of
% https://github.com/labgas/proj-emosymp/blob/main/prep/LaBGAS_extract_confound_reg_fMRIprep.m
% https://github.com/labgas/proj-emosymp/blob/main/firstlevel/LaBGAS_create_first_level_CANlab_tools_folder_structure.m
%
% It will extract noise regressors from fMRIprep output for
% the WAD-MRI study, including
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
% noise_regs & onsets files that can be loaded into CANlab DSGN structure
% (RECOMMENDED) using LaBGAS_get_firstlvl_dsgn_obj.m
% or directly into SPM first level batch (for the latter, use
% LaBGAS_first_level_batch_fMRIprep_conf.m)
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
% authors: Lukas Van Oudenhove & Iris Coppieters
% date:   April, 2021
%
%__________________________________________________________________________
% @(#)% WAD_prep_s3_extract_confound_reg_fMRIprep.m         v1.0        
% last modified: 2021/04/07


%% DEFINE AND CREATE DIRECTORIES, IDENTIFY SUBJECTS, AND SET OPTIONS
%--------------------------------------------------------------------------

% DIRECTORY STRUCTURE SETUP
% existing directories
rootdir = 'C:\Users\lukas\Dropbox (Dartmouth College)\4.Iris_WAD_MRI';
sourcedir = fullfile(rootdir,'sourcedata');
rawdir = fullfile(rootdir,'rawdata');
derivdir = fullfile(rootdir,'derivatives\fmriprep');

% create first level directory, and subdirs per model
% NOTE: these are not used further in this script, but by the next script,
% get_firstlvl_dsgn_obj.m
firstleveldir = fullfile(rootdir,'firstlevel');
if ~exist(firstleveldir,'dir')
    mkdir(firstleveldir);
end

models = {'model_1_pictures';'model_2_pictures_imagine';'model_3_pictures_fir'};
for i = 1:size(models,1)
    modeldir{i} = fullfile(firstleveldir,models{i});
    if ~exist(modeldir{i},'dir')
        mkdir(modeldir{i});
    end
end

% create subject list from derivdir, and check whether it matches subjdirs
% in rawdir, otherwise throw error
subjs = dir(fullfile(derivdir,'sub-*'));
idx = [subjs.isdir]';
subjs = {subjs(idx).name}';

subjs_raw = dir(fullfile(rawdir,'sub-*'));
idx_raw = [subjs_raw.isdir]';
subjs_raw = {subjs_raw(idx_raw).name}';

if ~isequal(subjs,subjs_raw)
    error('subject directories in derivatives and rawdata do not match; check and correct before proceeding');
end

clear idx idx_raw subjs_raw

% write subjectdirs in firstleveldir\modeldirs
sm=@(x)spm_mkdir(x); % defines spm_mkdir as an anonymous function sm
for j = 1:size(models,1)
    cd(modeldir{j});
    cellfun(sm,subjs);
end

% define runs
runs = {'run-A';'run-B'};

% SET OPTIONS
TR = 1;
spike_def = 'fMRIprep';
dvars_threshold = 2;
omit_spike_trials = 'no';
spike_additional_vols=0;
spikes_percent_threshold=0.15;

% THIS CHOICE OF OPTIONS CAN BE CONSIDERED LABGAS DEFAULTS, BUT MAY BE
% STUDY SPECIFIC, SO DISCUSS WITH LUKAS IF IN DOUBT!
% ALWAYS MAKE YOUR LOCAL COPY (OR BRANCH ON GITHUB) BEFORE MODIFYING THE CODE BELOW!


%% LOOP OVER SUBJECTS
%--------------------------------------------------------------------------

for sub=1:size(subjs,1)
    
    % DEFINE SUBJECT LEVEL DIRS
    subjrawdir=fullfile(rawdir,subjs{sub},'\ses-1\func');
    subjderivdir=fullfile(derivdir,subjs{sub},'\ses-1\func');
    cd(subjderivdir);
    cellfun(sm,runs);
    
    % DEFINE SUBJECT LEVEL FILENAMES
    % NOTE: Iris' raw images are unzipped .nii files, so no need to unzip
    % here
    rawimgs = dir(fullfile(subjrawdir,'*pic*bold.nii'));
    rawimgs = {rawimgs(:).name}';
    preprocimgs = dir(fullfile(subjderivdir,'s6-*pic*.nii.gz'));
    preprocimgs = {preprocimgs(:).name}';
    fmriprep_noisefiles = dir(fullfile(subjderivdir,'*task-pic*.tsv'));
    fmriprep_noisefiles = {fmriprep_noisefiles(:).name}';
    runnames = char(fmriprep_noisefiles);
    runnames = runnames(:,1:24); % this is study-specific, depends on length of subjectname and taskname - could be adapted based on regexp, but easy enough to adapt per study
        if ~isequal(size(rawimgs,1),size(preprocimgs,1),size(fmriprep_noisefiles,1)) 
            error('numbers of raw images, preprocessed images, and noise files do not match for %s, please check before proceeding',subjs{sub});
        end

    % CALCULATE AND EXTRACT CONFOUND REGRESSORS AND ONSETS, AND WRITE TO
    % FILE
    for run=1:size(fmriprep_noisefiles,1)
        
        % define subdir for this run
        rundir = fullfile(subjderivdir,runs{run});
        
        % move preprocessed image (unzip first) and fmriprep noisefile into rundir
        movefile(fullfile(subjderivdir,fmriprep_noisefiles{run}),fullfile(rundir,fmriprep_noisefiles{run}));
        gunzip(preprocimgs{run});
        niifile = dir(fullfile(subjderivdir,'s6-*pic*.nii'));
        movefile(fullfile(subjderivdir,niifile(1).name),fullfile(rundir,niifile(1).name));
        delete(preprocimgs{run});
        
        % CONFOUND REGRESSOR FILES
        % load confound regressor file generated by fMRIprep into Matlab table
        % variable
        R = readtable(fullfile(rundir,fmriprep_noisefiles{run}),'TreatAsEmpty','n/a','FileType', 'text', 'Delimiter', 'tab');

        % replace NaNs in first row with Os
        wh_replace = ismissing(R(1,:));
            if any(wh_replace)
                R{1, wh_replace} = zeros(1, sum(wh_replace)); % make array of zeros of the right size
            end

        % calculate and extract confound regressors
            if strcmpi(spike_def,'fMRIprep')==1 % switch would probably make more sense in this case, but this works too!

                % define regressors in fMRIprep output
                regs=R.Properties.VariableNames;
                spike_cols = contains(regs,'motion_outlier');
                Rspikes=R(:,spike_cols);
                Rspikes.spikes=sum(Rspikes{:,1:end},2);
                volume_idx = [1:height(R)]; 
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
                        additional_spikes_regs = zeros(height(R),size(spikes,2)*spike_additional_vols);
                            % This loop will create a separate column with ones in each row (TR) 
                            % we would like to consider a nuisance regressor
                            for i = 1:size(spikes,2) 
                                additional_spikes_regs(spikes(i)+1 : spikes(i)+spike_additional_vols,(i*spike_additional_vols-(spike_additional_vols-1)):(i*spike_additional_vols)) = eye(spike_additional_vols);
                            end
                        % if any spikes went beyond the end, trim it down
                        additional_spikes_regs = additional_spikes_regs(1:height(R),:);
                        % add the additional spikes to the larger matrix
                        R = [R array2table(additional_spikes_regs)];
                    end

                % remove redundant spike regressors
                regs = R.Properties.VariableNames;
                spike_cols = contains(regs,'motion_outlier');
                additional_spike_cols = contains(regs,'additional_spikes'); 
                [duplicate_rows, ~] = find(sum(R{:, spike_cols | additional_spike_cols}, 2)>1);
                    for i = 1:length(duplicate_rows) % This loop sets duplicate values to zero; drops them later (to keep indices the same during the loop)
                        [~,curr_cols] = find(R{duplicate_rows(i),:}==1);
                        R{duplicate_rows(i), curr_cols(2:end)} = 0;
                    end
                R = R(1:height(R), any(table2array(R)));

            elseif strcmpi(spike_def,'CANlab')==1

                % define raw image file
                raw_img_fname = rawimgs{run};

                % add in canlab spike detection (Mahalanobis distance)
                [g, mahal_spikes, gtrim, mahal_spikes_regs, snr] = scn_session_spike_id(fullfile(subjrawdir,raw_img_fname), 'doplot', 0); % CANlab function needs to be on your Matlab path
                delete('*.img'); % delete implicit mask .hdr/.img files generated by the CANlab function on the line above, since we don't need/use them
                delete('*.hdr');
                mahal_spikes_regs(:,1) = []; %drop gtrim which is the global signal
                R(:,contains(R.Properties.VariableNames,'motion_outlier'))=[]; % drop fmriprep motion outliers since we do not use them when spike_def = CANlab, and they cause redundancies
                R = [R array2table(mahal_spikes_regs)];

                % add in dvars spike regressors that are non-redundant with mahal spikes
                dvarsZ = [0; zscore(R.dvars(2:end))]; % first element of dvars always = 0, drop it from zscoring and set it to Z=0
                dvars_spikes = find(dvarsZ > dvars_threshold);
                same = ismember(dvars_spikes,mahal_spikes);
                dvars_spikes(same) = []; % drop the redundant ones
                dvars_spikes_regs = zeros(height(R),size(dvars_spikes,1));
                    for i=1:size(dvars_spikes,1)
                        dvars_spikes_regs(dvars_spikes(i),i) = 1;
                    end
                R = [R array2table(dvars_spikes_regs)];

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
                        additional_spikes_regs = additional_spikes_regs(1:height(R),:);
                        % add the additional spikes to the larger matrix
                        R = [R array2table(additional_spikes_regs)];
                    end

                % remove redundant spike regressors
                regs = R.Properties.VariableNames;
                spike_cols = contains(regs,'mahal_spikes'); 
                dvars_cols = contains(regs,'dvars_spikes'); 
                additional_spike_cols = contains(regs,'additional_spikes'); 

                [duplicate_rows, ~] = find(sum(R{:, spike_cols | dvars_cols | additional_spike_cols}, 2)>1);
                    for i = 1:size(duplicate_rows,1) %This loop sets duplicate values to zero; drops them later (to keep indices the same during the loop)
                        [~,curr_cols] = find(R{duplicate_rows(i),:}==1);
                        R{duplicate_rows(i), curr_cols(2:end)} = 0;
                    end
                R = R(1:size(mahal_spikes_regs,1), any(table2array(R)));
            else
                error('invalid spike_def option')
            end

        % Select confound and spike regressors to return for use in GLM 
        regs = R.Properties.VariableNames;
        motion_cols = contains(regs,'rot') | contains(regs,'trans');
        spike_cols = contains(regs,'mahal_spikes') | contains(regs,'motion_outlier'); 
        dvars_cols = contains(regs,'dvars_spikes'); 
        additional_spike_cols = contains(regs,'additional_spikes'); 
        Rselected = R(:,motion_cols | spike_cols | dvars_cols | additional_spike_cols);
        Rselected.csf = R.csf;
        Rspikes=R(:,spike_cols | dvars_cols | additional_spike_cols);
        Rspikes.spikes=sum(Rspikes{:,1:end},2);
        volume_idx = [1:height(R)]; 
        spikes = volume_idx(Rspikes.spikes==1)';

        % compute and output how many spikes total
        n_spike_regs = sum(dvars_cols | spike_cols | additional_spike_cols);
        n_spike_regs_percent = n_spike_regs / height(R);

        % print warning if #volumes identified as spikes exceeds
        % user-defined threshold
            if n_spike_regs_percent > spikes_percent_threshold
                warning('number of volumes identified as spikes exceeds threshold in %s',runnames(run,:))
            end

        % save confound regressors as .txt file
        filename_noise = strcat(rundir,'\noise_regs_',runnames(run,:));
        writetable(Rselected,filename_noise,'FileType','text','Delimiter','tab','WriteVariableNames',0);

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
        
    end % for loop runs
    
end % for loop subjects 