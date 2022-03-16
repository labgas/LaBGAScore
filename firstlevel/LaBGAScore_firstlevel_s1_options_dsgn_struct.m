%% LaBGAScore_firstlevel_s1_options_dsgn_struct.m
%
% This script sets the options and creates a CANlab-style DSGN structure
% variable, which are used by the subsequent script in the standard LaBGAS
% workflow, LaBGAScore_firstlevel_s2_fit_model, to fit first level models
% using CANlab and SPM functions
%
% USAGE
%
% Script should be run from the root directory of the superdataset, e.g.
% /data/proj_discoverie
% Script is highly study-specific, as conditions, contrasts, etc differ
% with design
% Hence, it is provided in LaBGAScore as an example and needs to be
% downloaded and adapted to the code subdataset for your study/project
% This example is from LaBGAS proj_erythritol_4a
% (https://gin.g-node.org/labgas/proj_erythritol_4a)
% 
%
% DEPENDENCIES
%
% LaBGAScore Github repo on Matlab path, with subfolders
% https://github.com/labgas/LaBGAScore
%
%
% OPTIONS
%
% MANDATORY OPTIONS
% 
% 1. spike_def
%
%   1. 'fMRIprep' uses spike regressors based on a combination of DVARS and FD thresholds 
%       set in fMRIprep arguments --fd-spike-threshold and --dvars-spike-threshold
%
%   2. 'CANlab' uses spike regressors based on CANlab's spike detection algorithm, which is based on Mahalanobis distance, and DVARS
%       function called: https://github.com/canlab/CanlabCore/blob/master/CanlabCore/diagnostics/scn_session_spike_id.m
%       for more info: https://github.com/canlab/CanlabScripts/blob/master/Scripts/Preprocessing/make_nuisance_covs_from_fmriprep_output2020.m
% 
% 2. omit_spike_trials
%
%   1. 'no' do not remove onsets of pain trials coinciding with a spike
%   2. 'yes' do remove
%       THIS IS NOT RECOMMENDED, WE PREFER TO DO THIS LATER BASED ON VIFS IN SINGLE TRIAL FIRST LEVEL ANALYSIS
%
% 3. spikes_percent_threshold
%
% set the maximum number of spikes (% of total volumes expressed as 0-1) you want to
% tolerate
%
% 4. vif_thresh
%
% set the threshold for variance inflation factor to flag regressors during
% model diagnosis
% default = 3, sensible range 1.3 (stringent) to 5 (lenient)
%
%
% NEEDED IF SPIKE_DEF = CANLAB
%
% 1. dvars_threshold
%
%   set the threshold for standardized dvars to define a spike
%   only used if spike_def = CANlab, otherwise set in fmriprep
%   --dvars-spike-threshold command
%   CANlab default is 3, but this is rather lenient, LaBGAS default 2
%
% 2. spike_additional_vols
%
%   set how many volumes after the spike you want to additionally regress out
%   be careful for task-based data since this quite aggressive approach is
%   mostly based on rs-fMRI, and beware of omitting too many volumes as well
%   as creating missingness not at random - THIS IS NOT RECOMMENDED
%
%
% NEEDED IF YOUR DESIGN INCLUDES PARAMETRIC MODULATORS
%
% 1. pmod_polynom
%
%   default = 1
%   can be changed to 2 or 3 for quadratic or cubic terms
%   be careful re orthogonalization issues in this case (see also below)
%   
% 2. pmod_name
%
%   variable name of the parametric modulator in events.tsv files
%
% 3. pmod_ortho_off
%
%   default = false
%   can be changed to true, but this is generally not recommended unless
%   you mean centered your pmods
%   for more info: 
%   https://andysbrainbook.readthedocs.io/en/latest/PM/PM_Overview.html
%   https://www.youtube.com/watch?v=iqwZmGOKqfM
%   CANlabReposGuide_Hackpad.pdf (see below)
%
%
% DSGN STRUCTURE
%
% The example below is maximally annotated, but more info on how to set up
% the DSGN structure and first level analyses with CANlab tools can be found at
%
% 1. canlab_glm_subject_levels('README') in Matlab terminal
% 2. canlab_glm_subject_levels('dsgninfo') in Matlab terminal
% 3. CANlabReposGuide_Hackpad.pdf
%       LaBGAS url: https://drive.google.com/drive/folders/1-M5UvibmsWXVCIrR31-qJNu506pDA_0t
%       CANlab url: https://drive.google.com/drive/folders/1G-_aDsylwOagCrS3ZMPPOmGsc2_eC4Nr
%
%__________________________________________________________________________
%
% authors: Lukas Van Oudenhove
% date:   March, 2022
%
%__________________________________________________________________________
% @(#)% LaBGAScore_firstlevel_s1_options_dsgn_struct.m         v1.0        
% last modified: 2022/03/16


%% SET OPTIONS
%--------------------------------------------------------------------------

% MANDATORY
spike_def = 'fMRIprep';
omit_spike_trials = 'no';
spikes_percent_threshold=0.15;
vif_thresh=3;

% ONLY NEEDED IF SPIKE_DEF = CANlab
dvars_threshold = 2;
spike_additional_vols=0;

% ONLY NEEDED IF YOU HAVE PARAMETRIC MODULATORS
pmod_polynom = 1; % polynomial expansion for pmods
pmod_name = 'rating'; % variable name of your pmod in events.tsv file
pmod_ortho_off = false; % turn off orthogonalization of pmods

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
    
    % flexible definition of contrasts - for more info
    % - help canlab_spm_contrast_job_luka
    % - https://github.com/canlab/CanlabCore/blob/master/CanlabCore/GLM_Batch_tools/canlab_glm_example_DSGN_setup.txt
    
    DSGN.regmatching = 'regexp'; % regular experession mode to match keywords you provide in cell arrays below with beta regressor names stored in the SPM.Vbeta.descrip field of your first level SPM.mat file
    % DSGN.defaultsuffix = '\*bf\(1\)$'; % adds this suffix to each keyword
    % DSGN.noscale = true; % default: false; if true, turn off scaling of contrasts to sum up to zero or one, not recommended if you have unequal amounts of run between subjects, unless you concatenate
    
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
    
    