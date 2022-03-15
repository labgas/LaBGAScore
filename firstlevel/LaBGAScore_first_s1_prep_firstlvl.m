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
% last modified: 2022/03/14
%
%
%% SET OPTIONS
%--------------------------------------------------------------------------

% MANDATORY
spike_def = 'fMRIprep';
omit_spike_trials = 'no';
spikes_percent_threshold=0.15;

% ONLY NEEDED IF SPIKE_DEF = CANlab
dvars_threshold = 2;
spike_additional_vols=0;

% ONLY NEEDED IF YOU HAVE PARAMETRIC MODULATORS
pmod_polynom = 1; % polynomial expansion for pmods
pmod_name = 'rating'; % variable name of your pmod in events.tsv file

% THIS CHOICE OF OPTIONS CAN BE CONSIDERED LABGAS DEFAULTS, BUT MAY BE
% STUDY SPECIFIC, SO DISCUSS WITH LUKAS IF IN DOUBT!


%% DEFINE DIRECTORIES AND RUNDIRNAMES
%--------------------------------------------------------------------------

% load standard BIDS directory structure from root dir

LaBGAScore_prep_s0_define_directories;

% define run directory names
% names are defaults, only change number

rundirnames = {'run-1';'run-2';'run-3';'run-4';'run-5';'run-6'};

% define rootdir for Github repos

githubrootdir = '/data/master_github_repos';


%% CREATE CANLAB DSGN STRUCTURE
%--------------------------------------------------------------------------    
% INPUT
    
    % REQUIRED FIELDS
    DSGN.metadata = "proj-erythritol_4a first level analysis model 1, i.e. modeling 4 conditions for sucrose, erythritol, sucralose, and water as long events (= duration of solution in mouth), with sweetness liking ratings as parametric modulators"; % field for annotation with study info, or whatever you like
    DSGN.modeldir = '/data/test_scripts/firstlevel/model_1_conds_pmods'; % directory where you want to write first level results for this model
    DSGN.subjects = derivsubjdirs';
    DSGN.funcnames = {'/func/run-1/s6*.nii',...
        '/func/run-2/s6*.nii',...
        '/func/run-3/s6*.nii',...
        '/func/run-4/s6*.nii',...
        '/func/run-5/s6*.nii',...
        '/func/run-6/s6*.nii'}; % cell array (one cell per session) of paths to functional files, relative to absolute path specific in DSGN.subjects
   
    % OPTIONAL FIELDS
    DSGN.concatenation = {}; % default: none; cell array of arrays of runs to concatenate; see documentation for when to concatenate, and how it works exactly
    DSGN.allowmissingfunc = true; % default: false; true will prevent erroring out when functional file is missing for at least one run is missing for at least one subject
    DSGN.customrunintercepts = {}; % default: none; will only work if DSGN.concatenation is specified; cell array of vectors specifying custom intercepts
    
% PARAMETERS

    DSGN.tr = 1.8; % repetition time (TR) in seconds
    DSGN.hpf = 180; % high pass filter in seconds; SPM default is 128, CANlab default is 180 since the brain response to pain stimuli last long and variance may be lost at shorter lengths, use scn_spm_design_check output, and review the SPM.mat in spm for diagnostics; 
    % STUDY-SPECIFIC: in this study, we stick with the CANlab default
    DSGN.fmri_t = 30; % microtime resolution - t=number of slices since we did slice timing; spm (and CANlab) default 16 can be kept for multiband w/o slice timing; TO BE CHECKED SINCE WE HAVE MULTIBAND WITH SLICE TIMING
    DSGN.fmri_t0 = 15; % microtime onset - reference slice used in slice timing correction; spm (and CANlab) default 1 can be kept for multiband w/o slice timing
    
% MODELING

    % REQUIRED FIELDS
    
    % cell array (one cell per session (i.e. run in our case)) of cell
    % arrays (one cell per condition) of MAT-file names, in fixed order:
    % all conditions of interest first, conditions of no interest last
    c=0;
    c=c+1;DSGN.conditions{c}={'sucrose' 'erythritol' 'sucralose' 'water' 'rating' 'swallow_rinse'};
    c=c+1;DSGN.conditions{c}={'sucrose' 'erythritol' 'sucralose' 'water' 'rating' 'swallow_rinse'};
    c=c+1;DSGN.conditions{c}={'sucrose' 'erythritol' 'sucralose' 'water' 'rating' 'swallow_rinse'};
    c=c+1;DSGN.conditions{c}={'sucrose' 'erythritol' 'sucralose' 'water' 'rating' 'swallow_rinse'};
    c=c+1;DSGN.conditions{c}={'sucrose' 'erythritol' 'sucralose' 'water' 'rating' 'swallow_rinse'};
    c=c+1;DSGN.conditions{c}={'sucrose' 'erythritol' 'sucralose' 'water' 'rating' 'swallow_rinse'};
    
    % OPTIONAL FIELDS
    
    % cell array (one cell per session) of cell arrays (one cell per condition) of cell arrays (one cell per modulator) of MAT-file names; set to {{}} if you don't want parametric modulators
    c=0;
    c=c+1;DSGN.pmods{c}={'liking_sucrose' 'liking_erythritol' 'liking_sucralose' 'liking_water'};
    c=c+1;DSGN.pmods{c}={'liking_sucrose' 'liking_erythritol' 'liking_sucralose' 'liking_water'};
    c=c+1;DSGN.pmods{c}={'liking_sucrose' 'liking_erythritol' 'liking_sucralose' 'liking_water'};
    c=c+1;DSGN.pmods{c}={'liking_sucrose' 'liking_erythritol' 'liking_sucralose' 'liking_water'};
    c=c+1;DSGN.pmods{c}={'liking_sucrose' 'liking_erythritol' 'liking_sucralose' 'liking_water'};
    c=c+1;DSGN.pmods{c}={'liking_sucrose' 'liking_erythritol' 'liking_sucralose' 'liking_water'};
%     DSGN.convolution; default hrf.derivs = [0 0]; structure specifying the convolution to use for conditions different fields required depending on convolution type; 
%     DSGN.ar1 = false; % autoregressive AR(1) to model serial correlations; SPM default is true, CANlab default is false, Tor recommends turning autocorrelation off, because this algorithm pools across the whole brain, and does not perform well in some situations; if you are performing a group analysis, the autocorrelation problem is not as concerning
    DSGN.notimemod = true; % default: false; if true, turn off time modulation of conditions, i.e. when you do not expect linear trends over time
%     DSGN.singletrials = {{}}; % a cell array (1 cell per session) of cell arrays (1 cell per condition) of (corresponding to DSGN.conditions) of true/false values indicating whether to convert specified condition to set of single trial conditions
%     DSGN.singletrialsall = false; % default: false; if true, set DSGN.singletrials to true for all conditions
    DSGN.modelingfilesdir = 'model_1_conds_pmods'; % name of subfolder which will be created within directory containing functional files where .mat files containing fields of DSGN structure will be saved
%     DSGN.allowemptycond = false; % default:false; if true, allow empty conditions
%     DSGN.allowmissingcondfiles = false; % default:false; if true, throw warning instead of error when no file(s) are found corresponding to a MAT-file name/wildcard
    DSGN.multireg = 'noise_regs'; % specify name for matfile with noise parameters you want to save
    
% CONTRASTS
    
    % OPTIONAL FIELDS
    
    % for flexible definition of contrasts - for more info
    % - help canlab_spm_contrast_job_luke
    % - https://github.com/canlab/CanlabCore/blob/master/CanlabCore/GLM_Batch_tools/canlab_glm_example_DSGN_setup.txt
    
    DSGN.regmatching = 'regexp'; % regular experession mode to match keywords you provide in cell arrays below with beta regressor names stored in the SPM.Vbeta.descrip field of your first level SPM.mat file
    % DSGN.defaultsuffix = '\*bf\(1\)$'; % adds this suffix to each keyword
    
    % REQUIRED FIELDS
    
    % cell array (one cell per contrast) of contrast definitions
    
    % unmodulated contrasts
    c=0;
    c=c+1;
    DSGN.contrasts{c} = {{'.*sucrose{1}\s[^x]'}}; % CON_0001; this regexp will select any beta regressor starting with "sucrose", followed by exactly one white space, but not followed by x - which is only the unmodulated regressors for the sucrose condition
    c=c+1;
    DSGN.contrasts{c} = {{'.*erythritol{1}\s[^x]'}}; % CON_0002
    c=c+1;
    DSGN.contrasts{c} = {{'.*sucralose{1}\s[^x]'}}; % CON_0003
    c=c+1;
    DSGN.contrasts{c} = {{'.*water{1}\s[^x]'}}; % CON_0004
    c=c+1;
    DSGN.contrasts{c} = {{'.*sucrose{1}\s[^x]'} {'.*water{1}\s[^x]'}}; % CON_0005
    c=c+1;
    DSGN.contrasts{c} = {{'.*erythritol{1}\s[^x]'} {'.*water{1}\s[^x]'}}; % CON_0006
    c=c+1;
    DSGN.contrasts{c} = {{'.*sucralose{1}\s[^x]'} {'.*water{1}\s[^x]'}}; % CON_0007
    c=c+1;
    DSGN.contrasts{c} = {{'.*sucrose{1}\s[^x]'} {'.*sucralose{1}\s[^x]'}}; % CON_0008
    c=c+1;
    DSGN.contrasts{c} = {{'.*sucrose{1}\s[^x]'} {'.*erythritol{1}\s[^x]'}}; % CON_0009
    c=c+1;
    DSGN.contrasts{c} = {{'.*erythritol{1}\s[^x]'} {'.*sucralose{1}\s[^x]'}}; % CON_0010
    
    % modulated contrasts
    c=c+1;
    DSGN.contrasts{c} = {{'.*liking_sucrose'}}; % CON_0011; this regexp will select any beta regressor with "liking_sucrose" anywhere in its name - which is only the modulated regressors for the sucrose condition
    c=c+1;
    DSGN.contrasts{c} = {{'.*liking_erythritol'}}; % CON_0012
    c=c+1;
    DSGN.contrasts{c} = {{'.*liking_sucralose'}}; % CON_0013
    c=c+1;
    DSGN.contrasts{c} = {{'.*liking_water'}}; % CON_0014
    c=c+1;
    DSGN.contrasts{c} = {{'.*liking_sucrose'} {'.*liking_water'}}; % CON_0015
    c=c+1;
    DSGN.contrasts{c} = {{'.*liking_erythritol'} {'.*liking_water'}}; % CON_0016
    c=c+1;
    DSGN.contrasts{c} = {{'.*liking_sucralose'} {'.*liking_water'}}; % CON_0017
    c=c+1;
    DSGN.contrasts{c} = {{'.*liking_sucrose'} {'.*liking_sucralose'}}; % CON_0018
    c=c+1;
    DSGN.contrasts{c} = {{'.*liking_sucrose'} {'.*liking_erythritol'}}; % CON_0019
    c=c+1;
    DSGN.contrasts{c} = {{'.*liking_erythritol'} {'.*liking_sucralose'}}; % CON_0020
    
    % OPTIONAL FIELDS
    
    % to define custom contrast names and weights
    % not needed strictly in this case, because this will be automatically generated for standard contrasts like this
    
    % unmodulated
    c=0;
    c=c+1;
    DSGN.contrastnames{c} = 'sucrose unmodulated'; % CON_0001
    DSGN.contrastweights{c} = [1];
    c=c+1;
    DSGN.contrastnames{c} = 'erythritol unmodulated'; % CON_0002
    DSGN.contrastweights{c} = [1];
    c=c+1;
    DSGN.contrastnames{c} = 'sucralose unmodulated'; % CON_0003
    DSGN.contrastweights{c} = [1];
    c=c+1;
    DSGN.contrastnames{c} = 'water unmodulated'; % CON_0004
    DSGN.contrastweights{c} = [1];
    c=c+1;
    DSGN.contrastnames{c} = 'sucrose unmodulated vs water unmodulated'; % CON_0005
    DSGN.contrastweights{c} = [1 -1];
    c=c+1;
    DSGN.contrastnames{c} = 'erythritol unmodulated vs water unmodulated'; % CON_0006
    DSGN.contrastweights{c} = [1 -1];
    c=c+1;
    DSGN.contrastnames{c} = 'sucralose unmodulated vs water unmodulated'; % CON_0007
    DSGN.contrastweights{c} = [1 -1];
    c=c+1;
    DSGN.contrastnames{c} = 'sucrose unmodulated vs sucralose unmodulated'; % CON_0008
    DSGN.contrastweights{c} = [1 -1];
    c=c+1;
    DSGN.contrastnames{c} = 'sucrose unmodulated vs erythritol unmodulated'; % CON_0009
    DSGN.contrastweights{c} = [1 -1];
    c=c+1;
    DSGN.contrastnames{c} = 'erythritol unmodulated vs sucralose unmodulated'; % CON_0010
    DSGN.contrastweights{c} = [1 -1];
    
    % modulated
    c=c+1;
    DSGN.contrastnames{c} = 'sucrose modulated'; % CON_0011
    DSGN.contrastweights{c} = [1];
    c=c+1;
    DSGN.contrastnames{c} = 'erythritol modulated'; % CON_0012
    DSGN.contrastweights{c} = [1];
    c=c+1;
    DSGN.contrastnames{c} = 'sucralose modulated'; % CON_0013
    DSGN.contrastweights{c} = [1];
    c=c+1;
    DSGN.contrastnames{c} = 'water modulated'; % CON_0014
    DSGN.contrastweights{c} = [1];
    c=c+1;
    DSGN.contrastnames{c} = 'sucrose modulated vs water modulated'; % CON_0015
    DSGN.contrastweights{c} = [1 -1];
    c=c+1;
    DSGN.contrastnames{c} = 'erythritol modulated vs water modulated'; % CON_0016
    DSGN.contrastweights{c} = [1 -1];    
    c=c+1;
    DSGN.contrastnames{c} = 'sucralose modulated vs water modulated'; % CON_0017
    DSGN.contrastweights{c} = [1 -1];
    c=c+1;
    DSGN.contrastnames{c} = 'sucrose modulated vs sucralose modulated'; % CON_0018
    DSGN.contrastweights{c} = [1 -1];
    c=c+1;
    DSGN.contrastnames{c} = 'sucrose modulated vs erythritol modulated'; % CON_0019
    DSGN.contrastweights{c} = [1 -1];
    c=c+1;
    DSGN.contrastnames{c} = 'erythritol modulated vs sucralose modulated'; % CON_0020
    DSGN.contrastweights{c} = [1 -1];
    

%--------------------------------------------------------------------------
% END OF STUDY-SPECIFIC CODE, BELOW SHOULD WORK OUT OF THE BOX
%--------------------------------------------------------------------------

%% MAKE SURE DEPENDENCIES ARE ON MATLAB PATH

% check whether necessary spm subdirs are on path, and add if needed

spmcanonicaldir = fullfile(spmrootdir,'canonical');
    if sum(contains(matlabpath,spmcanonicaldir)) == 0
        addpath(spmcanonicaldir,'-end');
        warning('adding %s to end of Matlab path',spmcanonicaldir)
    end
spmconfigdir = fullfile(spmrootdir,'config');
    if sum(contains(matlabpath,spmconfigdir)) == 0
        addpath(spmconfigdir,'-end');
        warning('adding %s to end of Matlab path',spmconfigdir)
    end
spmmatlabbatchdir = fullfile(spmrootdir,'matlabbatch');
    if sum(contains(matlabpath,spmmatlabbatchdir)) == 0
        addpath(spmmatlabbatchdir,'-end');
        warning('adding %s to end of Matlab path',spmmatlabbatchdir)
    end
spmtoolboxdir = fullfile(spmrootdir,'toolbox');
    if sum(contains(matlabpath,spmtoolboxdir)) == 0
        addpath(spmtoolboxdir,'-end');
        warning('adding %s to end of Matlab path',spmtoolboxdir)
    end
    
% check whether necessary CANlab Github repos are on Matlab path

  % CANLABCORE
    canlabcoredir = fullfile(githubrootdir,'CanlabCore');
        if ~isfolder(canlabcoredir) % canlabcore not yet cloned
          canlabcoreurl = "https://github.com/canlab/CanlabCore";
          canlabcoreclonecmd = ['git clone ' canlabcoreurl];
          cd(githubrootdir);
          [status,cmdout] = system(canlabcoreclonecmd);
          disp(cmdout);
              if status == -0
                  addpath(genpath(canlabcoredir,'-end'));
                  warning('git succesfully cloned %s to %s and added repo to Matlab path\n',canlabcoreurl, canlabcoredir)
              else
                  error('cloning %s into %s failed, please try %s in linux terminal before proceeding, or use Gitkraken\n',canlabcoreurl,canlabcoredir,canlabcoreclonecmd)
              end
          cd(rootdir);
          clear status cmdout
        elseif ~exist('fmri_data.m','file') % canlabcore cloned but not yet on Matlab path
            addpath(genpath(canlabcoredir,'-end'));
        end
        
  % CANLABPRIVATE
    canlabprivdir = fullfile(githubrootdir,'CanlabPrivate');
        if ~isfolder(canlabprivdir) % canlabprivate not yet cloned
          canlabprivurl = "https://github.com/canlab/CanlabPrivate";
          canlabprivclonecmd = ['git clone ' canlabprivurl];
          cd(githubrootdir);
          [status,cmdout] = system(canlabprivclonecmd);
          disp(cmdout);
              if status == -0
                  addpath(genpath(canlabprivdir,'-end'));
                  warning('git succesfully cloned %s to %s and added repo to Matlab path\n',canlabprivurl, canlabprivdir)
              else
                  error('cloning %s into %s failed, please try %s in linux terminal before proceeding, or use Gitkraken\n',canlabprivurl,canlabprivdir,canlabprivclonecmd)
              end
          cd(rootdir);
          clear status cmdout
        elseif ~exist('power_calc.m','file') % canlabprivate cloned but not yet on Matlab path
            addpath(genpath(canlabprivdir,'-end'));
        end

    
%% CREATE DIRECTORY STRUCTURE

% define mkdir as an anonymous function that can be applied to cell arrays
% to crease directory structure

sm=@(x)mkdir(x);

% create first level directory

firstleveldir = fullfile(rootdir,'firstlevel');
    if ~exist(firstleveldir,'dir')
        mkdir(firstleveldir);
    end

% create firstmodeldir

firstmodeldir = DSGN.modeldir;
    if ~exist(firstmodeldir,'dir')
        mkdir(firstmodeldir);
    end

% write subjectdirs in firstmodeldir

    if ~contains(ls(firstmodeldir),'sub-') % checks whether there are already subject dirs in firstmodeldir
        cd (firstmodeldir);
        cellfun(sm,derivsubjs);
    end
    
cd(rootdir);

% create list of subjectdirs in firstmodel dir

firstlist = dir(fullfile(firstmodeldir,'sub-*'));
firstsubjs = cellstr(char(firstlist(:).name));

    for firstsub = 1:size(firstsubjs,1)
        firstsubjdirs{firstsub,1} = fullfile(firstlist(firstsub).folder,firstlist(firstsub).name);
    end


%% LOOP OVER SUBJECTS
%--------------------------------------------------------------------------

for sub=1:size(derivsubjs,1)
    
    %% DEFINE SUBJECT LEVEL DIRS & FILENAMES
    
    subjderivdir = fullfile(derivsubjdirs{sub},'func');
    subjBIDSdir = fullfile(BIDSsubjdirs{sub},'func');
    subjfirstdir = firstsubjdirs{sub};
    
    BIDSimgs = dir(fullfile(subjBIDSdir,'*bold.nii.gz'));
    BIDSimgs = {BIDSimgs(:).name}';
    BIDSidx = ~contains(BIDSimgs,'rest'); % omit resting state scan if it exists
    BIDSimgs = {BIDSimgs{BIDSidx}}';
    
    derivimgs = dir(fullfile(subjderivdir,'s6-*.nii.gz'));
    derivimgs = {derivimgs(:).name}';
    derividx = ~contains(derivimgs,'rest'); % omit resting state scan if it exists
    derivimgs = {derivimgs{derividx}}';
    
    fmriprep_noisefiles = dir(fullfile(subjderivdir,'*desc-confounds_timeseries.tsv'));
    fmriprep_noisefiles = {fmriprep_noisefiles(:).name}';
    noiseidx = ~contains(fmriprep_noisefiles,'rest'); % omit resting state scan if it exists
    fmriprep_noisefiles = {fmriprep_noisefiles{noiseidx}}';
    
    % read events.tsv files with onsets, durations, and trial type
    eventsfiles = dir(fullfile(subjBIDSdir,'*events.tsv'));
    eventsfiles = {eventsfiles(:).name}';
    
    for runname = 1:size(fmriprep_noisefiles,1)
        subjrunnames{runname} = strsplit(fmriprep_noisefiles{runname},'_desc');
        subjrunnames{runname} = subjrunnames{runname}{1};
        subjrundirnames{runname} = subjrunnames{runname}(end-4:end);
    end
    subjrunnames = subjrunnames';
    subjrundirnames = subjrundirnames';
        
    % create rundirs in subjderivdir if needed
    if ~isfolder(fullfile(subjderivdir,rundirnames{1}))
        cd(subjderivdir);
        cellfun(sm,rundirnames);
    end
    
    cd(rootdir);
    
    % sanity check #1: number of images & noise/event files
    if ~isequal(size(BIDSimgs,1),size(derivimgs,1),size(fmriprep_noisefiles,1),size(eventsfiles,1)) 
        error('\n numbers of raw images, preprocessed images, noise, and events files do not match for %s, please check BIDSimgs, derivimgs, fmriprep_noisefiles, and eventsfiles variables before proceeding',derivsubjs{sub});
    else
        warning('\n numbers of raw images, preprocessed images, noise, and events files match for %s, continuing',derivsubjs{sub});
    end

    %% LOOP OVER RUNS: CREATE NOISE REGRESSORS, ONSETS, DURATIONS, AND PARAMETRIC MODULATORS
    
    for run=1:size(fmriprep_noisefiles,1)
        
        % DEFINE SUBDIR FOR THIS RUN
        rundir = fullfile(subjderivdir,subjrundirnames{run});
        
        % MOVE FMRIPREP NOISEFILE AND SMOOTHED IMAGE INTO RUNDIR IF
            if ~contains(ls(rundir),fmriprep_noisefiles{run})
                copyfile(fullfile(subjderivdir,fmriprep_noisefiles{run}),fullfile(rundir,fmriprep_noisefiles{run}));
            end

            if ~contains(ls(rundir),derivimgs{run})
                copyfile(fullfile(subjderivdir,derivimgs{run}),fullfile(rundir,derivimgs{run}));
                gunzip(fullfile(rundir,derivimgs{run}));
                delete(fullfile(rundir,derivimgs{run}));
            end

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

                % unzip & define raw image file
                gunzip(BIDSimgs{run}); % raw images are needed when spike_def = CANlab, which calls a function that is incompatible with .nii.gz, hence we unzip
                uBIDSimg = dir(fullfile(subjBIDSdir,'*bold.nii'));
                uBIDSimg = {uBIDSimg(:).name}';

                % add in canlab spike detection (Mahalanobis distance)
                [g, mahal_spikes, gtrim, mahal_spikes_regs, snr] = scn_session_spike_id(fullfile(subjBIDSdir,uBIDSimg), 'doplot', 0); % CANlab function needs to be on your Matlab path
                delete('*.img'); % delete implicit mask .hdr/.img files generated by the CANlab function on the line above, since we don't need/use them
                delete('*.hdr');
                delete(uBIDSimg); % delete unzipped image since we don't need it anymore and it eats up space
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
                    % set duplicate values to zero; drops them later (to keep indices the same during the loop)
                    for i = 1:size(duplicate_rows,1) 
                        [~,curr_cols] = find(Rfull{duplicate_rows(i),:}==1);
                        Rfull{duplicate_rows(i), curr_cols(2:end)} = 0;
                    end
                Rfull = Rfull(1:size(mahal_spikes_regs,1), any(table2array(Rfull)));
            else
                error('\n invalid spike_def option')
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
                warning('\n number of volumes identified as spikes exceeds threshold in %s',subjrunnames{run})
            end

        % save confound regressors as matrix named R for use in
        % SPM/CANlab GLM model tools
        R=table2array(R);

        % define and create subdir for model
        runmodeldir = fullfile(rundir,DSGN.modelingfilesdir);
            if ~exist(runmodeldir,'dir')
                mkdir(runmodeldir);
            end

        % write confound regressors
        filename_noise_regs = fullfile(runmodeldir,DSGN.multireg);
        save(filename_noise_regs,'R');

        clear R* filename_events

        % EVENTS FILES
        % read events.tsv files
        O = readtable(fullfile(subjBIDSdir,eventsfiles{run}),'FileType', 'text', 'Delimiter', 'tab');
        O.trial_type = categorical(O.trial_type);
        O.Properties.VariableNames(categorical(O.Properties.VariableNames) == pmod_name) = {'pmod'};
        
        % omit trials that coincide with spikes if that option is chosen
            if strcmpi(omit_spike_trials,'yes')==1
                same=ismember(O.onset,spikes); % identify trials for which onset coincides with spike
                O(same,:)=[]; % get rid of trials coinciding with spikes
            elseif strcmpi(omit_spike_trials,'no')==1
            else
                error('\n invalid omit_spike_trials option')
            end

        % sanity check #2: conditions
        cat_conds = reordercats(categorical(DSGN.conditions{run}));
        cat_conds = categories(cat_conds);
        cat_trial_type = cellstr(unique(O.trial_type));

            if ~isequal(cat_trial_type,cat_conds)
                error('\n conditions in DSGN structure do not match conditions in %s, please check before proceeding',fmriprep_noisefiles{run})
            else 
                warning('\n conditions in DSGN structure match conditions in %s, continuing',fmriprep_noisefiles{run})
            end

        % initialize structures for conditions
            for cond = 1:size(DSGN.conditions{1},2)
                cond_struct{cond} = struct('name',{DSGN.conditions{run}(cond)}, ...
                    'onset',{{[]}}, ...
                    'duration',{{[]}});
            end

            clear cond

        % fill structures with onsets and durations
            for trial = 1:size(O.trial_type,1)
                cond = 1;
                while cond < size(DSGN.conditions{run},2) + 1
                    switch O.trial_type(trial)
                        case DSGN.conditions{run}{cond}
                                cond_struct{cond}.onset{1} = [cond_struct{cond}.onset{1},O.onset(trial)];
                                cond_struct{cond}.duration{1} = [cond_struct{cond}.duration{1},O.duration(trial)];
                    end
                cond = cond + 1;
                end
                continue
            end

            clear cond

        % add pmods to structures for conditions of interest if specified in DSGN
            if isfield(DSGN,'pmods')
                for pmod = 1:size(DSGN.pmods{run},2)
                    cond_struct{pmod}.pmod = struct('name',{DSGN.pmods{run}(pmod)}, ...
                        'param',{{[]}}, ...
                        'poly',{{pmod_polynom}});
                end
                
                clear pmod
                
                for trial = 1:size(O.trial_type,1)
                    pmod = 1; 
                    while pmod < size(DSGN.pmods{run},2) + 1
                        switch O.trial_type(trial)
                            case DSGN.conditions{run}{pmod}
                               cond_struct{pmod}.pmod.param{1} = [cond_struct{pmod}.pmod.param{1},O.pmod(trial)];
                        end
                    pmod = pmod + 1;
                    end
                    continue
                end
                
                clear pmod
                
            end


        % sanity check #3: design info
        nii = dir(fullfile(rundir,'*.nii')).name;
        nii_hdr = read_hdr(fullfile(rundir,nii)); % reads Nifti header of smoothed image into a structure

            for cond = 1:size(DSGN.conditions{1},2)
                DesignTiming(1,cond) = (max(cond_struct{cond}.onset{1}) + cond_struct{cond}.duration{1}(1,end));
            end

            clear cond

        maxDesignTiming = max(DesignTiming);
        boldDuration = nii_hdr.tdim*DSGN.tr;

            if boldDuration < maxDesignTiming
                error('\n End of last condition (%s sec) exceeds BOLD duration (%s sec) in %s, please check before proceeding', num2str(maxDesignTiming), num2str(boldDuration), subjrunnames{run})
            else
                warning('\n End of last condition (%s sec) does not exceed BOLD duration (%s sec) in %s, continuing', num2str(maxDesignTiming), num2str(boldDuration), subjrunnames{run})
            end

        % save events file as .mat file
            for cond = 1:size(DSGN.conditions{1},2)
                struct = cond_struct{cond};
                save(fullfile(runmodeldir,char(cond_struct{cond}.name)),'-struct','struct');
                clear struct
            end
        
    end % for loop runs
    
    %% FIT FIRST LEVEL MODEL
    
    fprintf('\n Running on subject directory %s\n',DSGN.subjects{sub});
    canlab_glm_subject_levels(DSGN,'subjects',DSGN.subjects(sub),'overwrite','nolinks','noreview');
    
    
    %% DIAGNOSE FIRST LEVEL MODEL
    
    subjfirstdiagnosedir = fullfile(subjfirstdir,'diagnostics');
        if ~exist(subjfirstdiagnosedir,'dir')
            mkdir(subjfirstdiagnosedir);
        end
        
    cd(subjfirstdiagnosedir);
        
    diagnose_struct = struct('useNewFigure',false,'maxHeight',800,'maxWidth',1600,...
        'format','html','outputDir',subjfirstdiagnosedir,...
        'showCode',true);
    
    publish('LaBGAScore_first_s2_diagnose_firstlvl.m',diagnose_struct)
    delete('High_pass_filter_analysis.png','Variance_Inflation.png','LaBGAScore_first_s2_diagnose_firstlvl.png'); % getting rid of some redundant output images due to the use of publish()
    
    cd(rootdir);
    
    
end % for loop subjects 