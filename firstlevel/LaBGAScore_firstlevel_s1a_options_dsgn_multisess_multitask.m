%% LaBGAScore_firstlevel_s1a_options_dsgn_multisess_multitask.m
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
% LABGAS_OPTIONS STRUCTURE
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
%   set the maximum number of spikes (% of total volumes expressed as 0-1) you want to
%   tolerate
%
% 4. vif_thresh
%
%   set the threshold for variance inflation factor to flag regressors during
%   model diagnosis
%   default = 2, sensible range 1.3 (stringent) to 5 (lenient)
%
% 5. movement_reg_quadratic
%   
%   default = true
%   change to false to omit quadratic terms of movement parms and their
%   first-order derivatives
%
% 
% OPTIONAL
%
% 1. subjs2analyze
%
%   cell array of subjects in derivdir you want to analyze, empty cell array
%   if you want to loop over all subjects
%   NOTE: THIS OPTION IS NOT YET IMPLEMENTED IN THE NEXT SCRIPT, HENCE LEAVE
%           CELL ARRAY EMPTY OR COMMENT OUT FOR NOW
%
%
% SPIKE OPTIONS
%
% 1. dvars_threshold
%
%   set the threshold for standardized dvars to define a spike
%   MANDATORY if spike_def = CANlab, otherwise set in fmriprep
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
% MANDATORY ONLY IF YOUR DESIGN INCLUDES PARAMETRIC MODULATORS
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
%   should be changed to true if you have more than one pmod per condition,
%   which will cause the script to demean your pmods per condition, to
%   avoid the serial orthogonalization in spm
%   for more info: 
%   https://www.youtube.com/watch?v=iqwZmGOKqfM
%   CANlabReposGuide_Hackpad.pdf (see below)
%   https://www.bobspunt.com/resources/teaching/single-subject-analysis/parametric-modulation/#:~:text=In%20SPM%2C%20parametric%20modulators%20are,serial%20orthogonalisationto%20orthogonalise%20the%20parameters.
%   https://imaging.mrc-cbu.cam.ac.uk/imaging/ParametricModulations
%   https://www.researchgate.net/post/How_to_define_parametric_modulators_for_multiple_conditions_in_SPM
%
% 4. pmod_type
%   default = 'parametric_standard', alternative is 'parametric_singleregressor'
%   'standard' is for designs where each condition (of interest) has an
%   unmodulated and (orthogonalized/demeaned) modulated regressor
%   'singleregressor' is for designs where you only have modulated
%   regressors, which is less common, and beware of orthogonalization
%   issues
%   help(onsets2fmridesign) for more info about these options
%
%
% PLOTTING, AND THRESHOLDING & MASKING OPTIONS FOR DISPLAY PURPOSES
%  
% 1. plotdesign
%   default true, false suppresses design plots 
%   (saves time, but generally not recommended)
%
% 2. plotmontages
%   default true, false suppresses montages of first level con images
%   (saves a lot of time, but generally not recommended)
%   
%   inputs for the CANlab method statistic_image.threshold called by
%   LaBGAScore_firstlevel_s3_diagnose_model.m:
%
% 3. input_threshold
%   p-value or range of raw values, depending on thresh_type
%
% 4. thresh_type
%   threshold type which can be one of:
%     - 'fdr' : FDR-correct based on p-values already stored in image .p field
%     - 'bfr' : Bonferroni correction (FWE)
%     - 'unc' : Uncorrected p-value threshold: p-value, e.g., .05 or .001
%     - 'extent', 'cluster_extent' : Cluster extent correction with GRF at p < .05 corrected, primary threshold determined by input_threshold
%     - 'raw-between' : threshold raw image values; save those > input_threshold(1) and < input_threshold(2)
%     - 'raw-outside' : threshold raw image values; save those < input_threshold(1) or > input_threshold(2)
%
%   for more info:
%       doc statistic_image in Matlab terminal
%       edit statistic_image.threshold
%
% 5. k
%   extent threshold, in voxels
%
% 6. mask
%   image you want to mask with
% 
%
% DSGN STRUCTURE
%
%   The example below is maximally annotated, but more info on how to set up
%   the DSGN structure and first level analyses with CANlab tools can be found at
%
%   1. canlab_glm_subject_levels('README') in Matlab terminal
%   2. canlab_glm_subject_levels('dsgninfo') in Matlab terminal
%   3. CANlabReposGuide_Hackpad.pdf
%       LaBGAS url: https://drive.google.com/drive/folders/1-M5UvibmsWXVCIrR31-qJNu506pDA_0t
%       CANlab url: https://drive.google.com/drive/folders/1G-_aDsylwOagCrS3ZMPPOmGsc2_eC4Nr
%   
%   For info on contrast specification
%
%   1. help canlab_spm_contrast_job_luka in Matlab terminal
%   2. https://github.com/canlab/CanlabCore/blob/master/CanlabCore/GLM_Batch_tools/canlab_glm_example_DSGN_setup.txt
%
%__________________________________________________________________________
%
% authors: Lukas Van Oudenhove
% date:   September, 2023
%
%__________________________________________________________________________
% @(#)% LaBGAScore_firstlevel_s1a_options_dsgn_multisess_multitask.m    v1.1
% last modified: 2026/02/06
%
%
%% CREATE LABGAS_OPTIONS STRUCTURE
%--------------------------------------------------------------------------

% REQUIRED
LaBGAS_options.mandatory.spike_def = 'fMRIprep';
LaBGAS_options.mandatory.omit_spike_trials = 'no';
LaBGAS_options.mandatory.spikes_percent_threshold=0.15;
LaBGAS_options.mandatory.vif_thresh=2;
LaBGAS_options.movement_reg_quadratic = false; % change to false if you don't want to add quadratic terms for movement parameters and their first-order derivatives

% OPTIONAL
LaBGAS_options.subjs2analyze = {}; % enter subjects separated by comma if you only want to analyze selected subjects e.g. {'sub-01','sub-02'}; THIS IS NOT YET FULLY IMPLEMENTED HENCE LEAVE CELL ARRAY EMPTY OR COMMENT OUT OR DO NOT SPECIFY FIELD AT ALL

% SPIKE OPTIONS
LaBGAS_options.spikes.dvars_threshold = 2; % REQUIRED if spike_def = 'CANlab'
LaBGAS_options.spikes.spike_additional_vols=0; % OPTIONAL, NOT RECOMMENDED TO TURN ON

% OPTIONS FOR PLOTTTING, AND THRESHOLDING AND MASKING FIRST LEVEL IMAGES FOR DISPLAY
LaBGAS_options.display.plotdesign = true; % NOT RECOMMENDED TO TURN OFF
LaBGAS_options.display.plotmontages = true; % NOT RECOMMENDED TO TURN OFF
LaBGAS_options.display.input_threshold = 0.005;
LaBGAS_options.display.thresh_type = 'unc';
LaBGAS_options.display.k = 25;
LaBGAS_options.display.mask = which('gray_matter_mask_sparse.img');

% THE ABOVE OPTIONS CAN BE CONSIDERED LABGAS DEFAULTS, BUT MAY BE
% STUDY-SPECIFIC, SO DISCUSS WITH LUKAS IF IN DOUBT!

% REQUIRED IF YOU HAVE PARAMETRIC MODULATORS
LaBGAS_options.pmods.pmod_polynom = 1;
LaBGAS_options.pmods.pmod_name = 'wanting';
LaBGAS_options.pmods.pmod_ortho_off = false;
LaBGAS_options.pmods.pmod_type = 'parametric_standard';

if strcmpi(LaBGAS_options.mandatory.spike_def,'CANlab')==1 && ~isfield(LaBGAS_options.spikes,'dvars_threshold')
    error('spike_def option %s requires specification of LaBGAS_options.spikes.dvars_threshold, please specify before proceeding',LaBGAS_options.mandatory.spike_def)
end


%% DEFINE DIRECTORIES AND RUNDIRNAMES
%--------------------------------------------------------------------------

% load standard BIDS directory structure from root dir
% STUDY-SPECIFIC: replace LaBGAScore with name of study-specific script

LaBGAScore_prep_s0_define_directories;

% define run directory names
% STUDY-SPECIFIC: names are defaults, only change number
% NOTE: if you only have one task, use 'run-1', 'run-2', etc
%       if you have more than one task, use one model/script per task, and
%       add _<taskname> (as in your filenames in derivdir) to rundirnames

rundirnames = {'run-1_FID';'run-2_FID'};

% set number of sessions

nr_sess = 2;

% define tasknames (as in BIDS-compliant naming of your functional files in BIDSdir)

tasknames = {'food_images';'FID'};
taskname = tasknames{1}; % task for which you want to run firstlevel analysis

% define rootdir for Github repos
% MACHINE-SPECIFIC: change to your local path if not working on LaBGAS server

githubrootdir = '/data/master_github_repos';


%% CREATE CANLAB DSGN STRUCTURE
%--------------------------------------------------------------------------    
% INPUT
    
    % REQUIRED FIELDS
    DSGN.metadata = "proj-bitter-reward first level analysis model with FID task, i.e. modeling 6 conditions for cues with 0,2 or 10 win, and feedback with 0,2 or 10 win"; % field for annotation with study info, or whatever you like
    DSGN.modeldir = '/data/proj_bitter-reward/firstlevel/model_2_FID'; % directory where you want to write first level results for this model
        if ~isfield(LaBGAS_options,'subjs2analyze')
            DSGN.subjects = derivsubjdirs';
        elseif ~isempty(LaBGAS_options.subjs2analyze)
            [C,~,~] = intersect(derivsubjs,LaBGAS_options.mandatory.subjs2analyze);
            if ~isequal(C',LaBGAS_options.mandatory.subjs2analyze)
                error('\n subject %s defined in LaBGAS_options.mandatory.subjs2analyze not present in %s, please check before proceeding',LaBGAS_options.mandatory.subj2analyze{~ismember(LaBGAS_options.mandatory.subjs2analyze,C)},derivdir);
            else
                DSGN.subjects = cell(1,size(LaBGAS_options.mandatory.subjs2analyze,2));
                    for sub = 1:size(DSGN.subjects,2)
                        DSGN.subjects{sub} = fullfile(derivdir,LaBGAS_options.mandatory.subjs2analyze{sub});
                    end
            end
        else
            DSGN.subjects = derivsubjdirs';
        end
    DSGN.funcnames = {['/ses-1/func/' rundirnames{1} '/s6*ses-1*' taskname '*run-1*.nii'],...
        ['/ses-1/func/' rundirnames{2} '/s6*ses-1*' taskname '*run-2*.nii'],...
        ['/ses-2/func/' rundirnames{1} '/s6*ses-2*' taskname '*run-1*.nii'],...
        ['/ses-2/func/' rundirnames{2} '/s6*ses-2*' taskname '*run-2*.nii']}; % cell array (one cell per session) of paths to functional files, relative to absolute path specified in DSGN.subjects
   
    % OPTIONAL FIELDS
%     DSGN.concatenation = {[1:6]}; % default: none; cell array of arrays of runs to concatenate; see documentation for when to concatenate, and how it works exactly
    DSGN.allowmissingfunc = true; % default: false; true will prevent erroring out when functional file is missing for at least one run is missing for at least one subject
%     DSGN.customrunintercepts = {1:6}; % default: none; will only work if DSGN.concatenation is specified; cell array of vectors specifying custom intercepts, NOT YET FULLY TESTED 
    
% PARAMETERS

    DSGN.tr = 2.5; % repetition time (TR) in seconds
    DSGN.hpf = 180; % high pass filter in seconds; SPM default is 128, CANlab default is 180 since the brain response to pain stimuli last long and variance may be lost at shorter lengths, use scn_spm_design_check output, and review the SPM.mat in spm for diagnostics; 
    % STUDY-SPECIFIC: in this study, we stick with the CANlab default
    DSGN.fmri_t = 46; % microtime resolution - t=number of slices since we did slice timing; spm (and CANlab) default 16 can be kept for multiband w/o slice timing; TO BE CHECKED SINCE WE HAVE MULTIBAND WITH SLICE TIMING
    DSGN.fmri_t0 = 23; % microtime onset - reference slice used in slice timing correction; spm (and CANlab) default 1 can be kept for multiband w/o slice timing
    
% MODELING

    % REQUIRED FIELDS
    
    % cell array (one cell per session (i.e. run in our case)) of cell
    % arrays (one cell per condition) of MAT-file names, in fixed order:
    % all conditions of interest first, conditions of no interest last
    
    c=0;
    c=c+1;DSGN.conditions{c}={'bit_cue_C0','bit_cue_C2','bit_cue_C10','bit_feedback_C0','bit_feedback_C2_win','bit_feedback_C10_win'};
    c=c+1;DSGN.conditions{c}={'bit_cue_C0','bit_cue_C2','bit_cue_C10','bit_feedback_C0','bit_feedback_C2_win','bit_feedback_C10_win'};
    c=c+1;DSGN.conditions{c}={'pla_cue_C0','pla_cue_C2','pla_cue_C10','pla_feedback_C0','pla_feedback_C2_win','pla_feedback_C10_win' };
    c=c+1;DSGN.conditions{c}={'pla_cue_C0','pla_cue_C2','pla_cue_C10', 'pla_feedback_C0','pla_feedback_C2_win','pla_feedback_C10_win'};
    
    
    % OPTIONAL FIELDS
    
    % cell array (one cell per session) of cell arrays (one cell per condition) of cell arrays (one cell per modulator) of MAT-file names; set to {{}} if you don't want parametric modulators
%   c=0;
%   c=c+1;DSGN.conditions{c}={};
%   c=c+1;DSGN.conditions{c}={};
%   c=c+1;DSGN.conditions{c}={};
%   c=c+1;DSGN.conditions{c}={};
%   DSGN.convolution.type; default hrf, which means canonical hrf - other options: fir, spline (the latter is not yet implemented @LaBGAS, help needed from Tor/Martin/Bogdan)
%   DSGN.convolution.time; default 0, which means no time derivative
%   DSGN.convolution.dispersion: default 0, which means no dispersion derivative
%   DSGN.ar1 = false; % autoregressive AR(1) to model serial correlations; SPM default is true, CANlab default is false, Tor recommends turning autocorrelation off, because this algorithm pools across the whole brain, and does not perform well in some situations; if you are performing a group analysis, the autocorrelation problem is not as concerning
    DSGN.notimemod = true; % default: false; if true, turn off time modulation of conditions, i.e. when you do not expect linear trends over time
%     DSGN.singletrials = {{}}; % a cell array (1 cell per session) of cell arrays (1 cell per condition) of (corresponding to DSGN.conditions) of true/false values indicating whether to convert specified condition to set of single trial conditions
%     DSGN.singletrialsall = false; % default: false; if true, set DSGN.singletrials to true for all conditions
    DSGN.modelingfilesdir = 'model_2_FID'; % name of subfolder which will be created within directory containing functional files where .mat files containing fields of DSGN structure will be saved
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
    
    c=0;
    c=c+1;
    DSGN.contrasts{c} = {{'.*bit_cue_C0{1}\s[^x]'}}; % CON_0001; this regexp will select any beta regressor starting with "bit_cue", followed by exactly one white space, but not followed by x - which is only the unmodulated regressors for the cue_0 condition
    c=c+1;
    DSGN.contrasts{c} = {{'.*bit_cue_C2{1}\s[^x]'}}; % CON_0002
    c=c+1;
    DSGN.contrasts{c} = {{'.*bit_cue_C10{1}\s[^x]'}}; % CON_0003
    c=c+1;
    DSGN.contrasts{c} = {{'.*pla_cue_C0{1}\s[^x]'}}; % CON_0004
    c=c+1;
    DSGN.contrasts{c} = {{'.*pla_cue_C2{1}\s[^x]'}}; % CON_0005
    c=c+1;
    DSGN.contrasts{c} = {{'.*pla_cue_C10{1}\s[^x]'}}; % CON_0006
    c=c+1;
    DSGN.contrasts{c} = {{'.*bit_feedback_C0{1}\s[^x]'}}; % CON_007
    c=c+1;
    DSGN.contrasts{c} = {{'.*bit_feedback_C2_win{1}\s[^x]'}}; % CON_008
    c=c+1;
    DSGN.contrasts{c} = {{'.*bit_feedback_C10_win{1}\s[^x]'}}; % CON_009
    c=c+1;
    DSGN.contrasts{c} = {{'.*pla_feedback_C0{1}\s[^x]'}}; % CON_010 
    c=c+1;
    DSGN.contrasts{c} = {{'.*pla_feedback_C2_win{1}\s[^x]'}}; % CON_011
    c=c+1;
    DSGN.contrasts{c} = {{'.*pla_feedback_C10_win{1}\s[^x]'}}; % CON_012
    c=c+1;
    DSGN.contrasts{c} = {{'.*bit_cue_C10{1}\s[^x]'} {'.*bit_cue_C0{1}\s[^x]'}}; % CON_013;
    c=c+1;
    DSGN.contrasts{c} = {{'.*pla_cue_C10{1}\s[^x]'} {'.*pla_cue_C0{1}\s[^x]'}}; % CON_014;
    c=c+1;
    DSGN.contrasts{c} = {{'.*bit_feedback_C10_win{1}\s[^x]'} {'.*bit_feedback_C0{1}\s[^x]'}}; % CON_015;
    c=c+1;
    DSGN.contrasts{c} = {{'.*pla_feedback_C10_win{1}\s[^x]'} {'.*pla_feedback_C0{1}\s[^x]'}}; % CON_016;
    c=c+1;
    DSGN.contrasts{c} = {{'.*bit_cue_C10{1}\s[^x]'} {'.*bit_cue_C0{1}\s[^x]'} {'.*pla_cue_C10{1}\s[^x]'} {'.*pla_cue_C0{1}\s[^x]'}}; % CON_017;
    c=c+1;
    DSGN.contrasts{c} = {{'.*bit_feedback_C10_win{1}\s[^x]'} {'.*bit_feedback_C0{1}\s[^x]'} {'.*pla_feedback_C10_win{1}\s[^x]'} {'.*pla_feedback_C0{1}\s[^x]'}}; % CON_018;
    
   
    % OPTIONAL FIELDS
    
    % to define custom contrast names and weights
    % not needed strictly in this case, because this will be automatically generated for standard contrasts like this
    
    c=0;
    c=c+1;
    DSGN.contrastnames{c} = 'bitter cue 0 SP'; % CON_0001
    DSGN.contrastweights{c} = [1];
    c=c+1;
    DSGN.contrastnames{c} = 'bitter cue 2 SP'; % CON_0002
    DSGN.contrastweights{c} = [1];
    c=c+1;
    DSGN.contrastnames{c} = 'bitter cue 10 SP'; % CON_0003
    DSGN.contrastweights{c} = [1];
    c=c+1;
    DSGN.contrastnames{c} = 'placebo cue 0 SP'; % CON_0004
    DSGN.contrastweights{c} = [1];
    c=c+1;
    DSGN.contrastnames{c} = 'placebo cue 2 SP'; % CON_0005
    DSGN.contrastweights{c} = [1];
    c=c+1;
    DSGN.contrastnames{c} = 'placebo cue 10 SP'; % CON_0006
    DSGN.contrastweights{c} = [1];
    c=c+1;
    DSGN.contrastnames{c} = 'bitter feedback 0 SP'; % CON_0007
    DSGN.contrastweights{c} = [1];
    c=c+1;
    DSGN.contrastnames{c} = 'bitter feedback 2 SP'; % CON_0008
    DSGN.contrastweights{c} = [1];
    c=c+1;
    DSGN.contrastnames{c} = 'bitter feedback 10 SP'; % CON_0009
    DSGN.contrastweights{c} = [1];
    c=c+1;
    DSGN.contrastnames{c} = 'placebo feedback 0 SP'; % CON_0010
    DSGN.contrastweights{c} = [1];
    c=c+1;
    DSGN.contrastnames{c} = 'placebo feedback 2 SP'; % CON_0011
    DSGN.contrastweights{c} = [1];
    c=c+1;
    DSGN.contrastnames{c} = 'placebo feedback 10 SP'; % CON_0012
    DSGN.contrastweights{c} = [1];
    c=c+1;
    DSGN.contrastnames{c} = 'bitter contrast cue 10 SP minus cue 0 SP'; % CON_0013
    DSGN.contrastweights{c} = [1 -1];
    c=c+1;
    DSGN.contrastnames{c} = 'placebo contrast cue 10 SP minus cue 0 SP'; % CON_0014
    DSGN.contrastweights{c} = [1 -1]; 
    c=c+1;
    DSGN.contrastnames{c} = 'bitter feedback 10 SP minus feedback 0 SP'; % CON_0015
    DSGN.contrastweights{c} = [1 -1];
    c=c+1;
    DSGN.contrastnames{c} = 'placebo feedback 10 SP minus feedback 0 SP'; % CON_0016
    DSGN.contrastweights{c} = [1 -1];
    c=c+1;
    DSGN.contrastnames{c} = 'bitter vs placebo cue 10 SP vs cue 0 SP'; % CON_0017
    DSGN.contrastweights{c} = [1 -1 -1 1];
    c=c+1;
    DSGN.contrastnames{c} = 'bitter vs placebo feedback 10 SP minus feedback 0 SP'; % CON_0018
    DSGN.contrastweights{c} = [1 -1 -1 1];