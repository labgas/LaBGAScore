%% LaBGAScore_decoding_template_xclass_acc.m
%
%
% *USAGE*
%
% WITHIN-SUBJECT (CROSS-)CLASSIFICATION ACCURACY WITH LEAVE-ONE-RUN-OUT CV
%
% This script uses The Decoding Toolbox
% <https://sites.google.com/site/tdtdecodingtoolbox/ (TDT)> to classify among
% conditions WITHIN each subject, using that subject's first-level betas, and
% then summarizes the per-subject accuracies at group level.
%
% Note the level: the unit of classification is a condition beta within a
% subject, and the cross-validation unit is a RUN. This is a different question
% from LaBGAScore_decoding_SVM_between_subjects.m, which classifies SUBJECTS
% into groups; the two are not comparable measurements and their accuracies
% should not be reported side by side as if they were. See
% decoding_toolbox/README.md
%
% The script executes the following steps:
%
% 1. Prep work
%    - Calls <study_prefix>_prep_s0_define_directories, the first level
%      s1_options_dsgn_struct (for DSGN.modeldir) and the second level
%      a_set_up_paths_always_run_first (for resultsdir/htmlsavedir) if the
%      variables they define are not already in the workspace.
%    - Lists the fitted subjects under DSGN.modeldir and creates the TDT
%      output directories at first and second level.
%
% 2. Per subject, over the loop
%    - Loads that subject's SPM.mat.
%    - Resamples the mask to the first subject's beta space and writes it once,
%      then CHECKS every later subject against that reference.
%      NOTE: one resampled mask is only valid while every subject shares a
%          space; a mismatch is now an error rather than a silent misalignment
%    - Reads the condition names from SPM.Sess(1).U, by POSITION in
%      conds2include.
%    - Builds the leave-one-run-out design and runs the decoding.
%      NOTE: TDT functions called here
%          a) decoding_defaults      - cfg with TDT's defaults
%          b) design_from_spm        - beta names and run numbers from SPM.mat
%          c) decoding_describe_data - labels, chunks and file list into cfg
%          d) make_design_cv         - the leave-one-run-out CV design
%          e) decoding               - runs it
%    - Collects the confusion matrix and plots it.
%
% 3. Group level
%    - Mean, SD, median and IQR confusion matrices across subjects.
%    - Repeated-measures ANOVA across conditions (fitrm/ranova) with
%      Tukey-Kramer corrected post-hoc comparisons.
%    - Wilcoxon sign-rank test of each condition's accuracy against chance.
%      NOTE: the ANOVA is PARAMETRIC while the sign-rank test and the
%          median/IQR reporting are not. Both are reported deliberately; if
%          the accuracies are visibly non-normal, lead with the
%          non-parametric results
%      NOTE: chance level depends on the SCALE of TDT's confusion_matrix
%          output. The script derives it from the data (percent vs
%          proportion) and prints which one it used - check that line
%    - Saves results and figures under <resultsdir>/TDT/<analysis_name>.
%
% Run this script with Matlab's publish function to generate an html report:
% publish('LaBGAScore_decoding_template_xclass_acc.m','outputDir',htmlsavedir)
% PREFER LaBGAScore_prov_publish (clean/), which additionally records the
% commit of every dependency the run reached - see clean/README_provenance.md
% NOTE: publish() catches a script error INTO the html and returns normally,
%     so a crashed run raises nothing and exits 0. clean/LaBGAScore_run_reports.m
%     exists to close that hole
%
%
% *OPTIONS*
%
% Set in the "SET OPTIONS" section below:
%
% * results_suffix                  string to add to results files, useful if you want to run multiple analyses for example with different masks or conditions
%
% * conds2include                   indices of conditions in DSGN.conditions (or SPM.Sess(x).U.name); assumes all conditions are present in all runs!
%
% * mask_name                       absolute path to mask, or name of mask file if already on Matlab path
%
% * scaling_regime                  'kernel' (TDT's recommended estimation = 'all', label-blind and faster) or
%                                   'strict' (train-fold-only 'across' scaling). Same choice, same reasoning,
%                                   as in LaBGAScore_decoding_SVM_between_subjects.m
%
% Taken from the workspace rather than set here, so the scripts that define
% them must have run (this script calls each if its variable is missing):
%
% * DSGN.modeldir                   first level model dir, from <study_prefix>_firstlevel_<mx>_s1_options_dsgn_struct
% * rootdir, githubrootdir          from <study_prefix>_prep_s0_define_directories
% * resultsdir, htmlsavedir         from a_set_up_paths_always_run_first
%
%
% *DEPENDENCIES*
%
% 1. The Decoding Toolbox (TDT) on your Matlab path
%   https://sites.google.com/site/tdtdecodingtoolbox/
% 2. spm12 on your Matlab path - first level SPM.mat and betas are the input
% 3. CANlab's CanlabCore on your Matlab path - fmri_mask_image, fmri_data,
%   resample_space
% 4. LaBGAScore Github repo on your Matlab path, with subfolders
%   https://github.com/labgas/LaBGAScore
% 5. Matlab's Statistics and Machine Learning Toolbox - fitrm, ranova,
%   multcompare, signrank
%
% For more info on TDT, check out the
% <https://www.frontiersin.org/articles/10.3389/fninf.2014.00088/full
% accompanying paper> and/or
% <https://andysbrainbook.readthedocs.io/en/latest/ML/ML_Short_Course/ML_05_Haxby_MVPA.html
% Andy Jahn's tutorials>
%
%
% *NOTES*
%
% INPUTS: first level betas and SPM.mat per subject under DSGN.modeldir; a mask
%
% OUTPUTS: per-subject TDT results under <subject>/TDT/<analysis_name>/, and
% group level confusion matrices, stats and figures under
% <resultsdir>/TDT/<analysis_name>/
%
% ASSUMPTIONS, none of which are checked beyond what is noted above:
%
% 1. ALL CONDITIONS ARE PRESENT IN ALL RUNS. conds2include indexes
%    SPM.Sess(1).U, and the same condition structure is assumed for every run
%    and every subject.
% 2. Condition NAMES are read from session 1 only. In a multi-session design
%    where the sessions differ, this is wrong.
% 3. Every subject is in the same space as subject 1 (now checked, and an
%    error if not).
% 4. SPM.Sess(1).U(k).name carries a trailing blank space relative to the
%    names in DSGN.conditions - hence the char() conversion.
%
% This script was adapted from TDT's own decoding_template.m, which is why the
% per-subject section still carries that template's commentary.
%
% -------------------------------------------------------------------------
%
% modified by: Lukas Van Oudenhove (adapted from decoding_template.m)
%
% date:   KU Leuven, June, 2023
%
% -------------------------------------------------------------------------
%
% LaBGAScore_decoding_template_xclass_acc.m         v2.0
%
% last modified: 2026/09/30
%
%
%
%% SET OPTIONS
% -------------------------------------------------------------------------

results_suffix = 'all_conds';
conds2include = 1:4; % indices of conditions in DSGN.conditions (or SPM.Sess(x).U.name); assumes all conditions are present in all runs!
mask_name = 'gm_mask_canlab2023_coarse_fmriprep20_0_20.nii'; % absolute path to mask, or name of mask file if already on Matlab path

% scaling regime, as in LaBGAScore_decoding_SVM_between_subjects.m. This is a
% methodological choice, not a speed knob:
%   'kernel' - TDT's own recommendation, scaling estimated on ALL data. Label-
%              blind, so it cannot carry condition information across the
%              train/test split, and it is much faster.
%   'strict' - scaling estimated on the TRAINING FOLD ONLY ('across'). Slower,
%              and the conservative choice if you would rather not estimate any
%              quantity on data the classifier is about to be tested on.
scaling_regime = 'kernel'; % 'kernel' | 'strict'


%% PREP WORK
% -------------------------------------------------------------------------

% set analysis name

analysis_name = 'xclass_acc';

% check whether LaBGAScore_prep_s0_define_directories has been run
% STUDY-SPECIFIC: replace LaBGAScore with study name in code below

if ~exist('rootdir','var') || ~exist('githubrootdir','var')
    warning('\nrootdir and/or githubrootdir variable not found in Matlab workspace, running LaBGAScore_prep_s0_define_directories before proceeding')
    LaBGAScore_prep_s0_define_directories;
    cd(rootdir);
else
    cd(rootdir);
end

% check whether LaBGAScore_firstlevel_s1_options_dsgn_struct.m has been run
% STUDY-SPECIFIC: replace LaBGAScore with study name and add model index in code below

if ~exist('DSGN','var')
    warning('\nDSGN variable not found in Matlab workspace, running LaBGAScore_firstlevel_s1_options_dsgn_struct.m before proceeding')
    LaBGAScore_firstlevel_s1_options_dsgn_struct;
end

% check whether LaBGAScore_firstlevel_s0_a_set_up_paths_always_run_first.m has been run
% STUDY-SPECIFIC: replace LaBGAScore with study name and add model index in code below

if ~exist('htmlsavedir','var')
    warning('\nhtmlsavedir variable not found in Matlab workspace, running a_set_up_paths_always_run_first before proceeding')
    a_set_up_paths_always_run_first;
end

% get first level dir info

firstmodeldir = DSGN.modeldir;

firstlist = dir(fullfile(firstmodeldir,'sub-*'));
firstlist = firstlist([firstlist(:).isdir]); % directories only, so a stray sub-* file is not taken for a subject
firstsubjs = cellstr(char(firstlist(:).name));

if isempty(firstlist)
    error('\nno sub-* directories found in %s, please check that the first level model has been fitted before proceeding',firstmodeldir);
end

% nsubs is used everywhere below. It used to be the LEAKED COUNTER of the loop
% that follows (firstsub), which happened to hold the right value but was
% undefined whenever that loop did not run.
nsubs = size(firstsubjs,1);
firstsubjdirs = cell(nsubs,1);

    for k = 1:nsubs
        firstsubjdirs{k,1} = fullfile(firstlist(k).folder,firstlist(k).name);
    end
clear k
    
clear firstlist

% set firstlevel TDT dir info & create dirs if needed

firstmodelTDTdir = fullfile(firstmodeldir,'TDT');
if ~isfolder(firstmodelTDTdir)
    mkdir(firstmodelTDTdir);
end

firstmodelTDTdir_mask = fullfile(firstmodelTDTdir,'mask');
if ~isfolder(firstmodelTDTdir_mask)
    mkdir(firstmodelTDTdir_mask);
end

% set secondlevel TDT dir info & create dir if needed

secondmodelTDTdir = fullfile(resultsdir,'TDT');
if ~isfolder(secondmodelTDTdir)
    mkdir(secondmodelTDTdir);
end

secondmodelTDTanalysisdir = fullfile(secondmodelTDTdir, analysis_name);
if ~isfolder(secondmodelTDTanalysisdir)
    mkdir(secondmodelTDTanalysisdir);
end

% pre-allocate loop vars

spmdotmats = cell(1,nsubs);
results_combined = cell(1,nsubs);
group_results = cell(1,size(conds2include,2));

for cond = 1:size(group_results,2)
    group_results{cond} = [];
end

confusion_matrices = [];


%% LOOP TDT TEMPLATE CODE OVER SUBJECTS
% -------------------------------------------------------------------------

for sub = 1:nsubs
    
    fprintf('\n\n');
    printhdr(['SUBJECT #', num2str(sub)]);
    fprintf('\n\n');

    spmdotmats{sub} = load(fullfile(firstsubjdirs{sub},'SPM.mat'));

    % This script is a template that can be used for a decoding analysis on 
    % brain image data. It is for people who have betas available from an 
    % SPM.mat and want to automatically extract the relevant images used for
    % classification, as well as corresponding labels and decoding chunk numbers
    % (e.g. run numbers). If you don't have this available, then use
    % decoding_template_nobetas.m

    % Make sure the decoding toolbox and your favorite software (SPM or AFNI)
    % are on the Matlab path (e.g. addpath('/home/decoding_toolbox') )
    % TDT
    % addpath('$ADD FULL PATH TO TDT TOOLBOX AS STRING OR MAKE THIS LINE A COMMENT IF IT IS ALREADY$')
    % assert(~isempty(which('decoding_defaults.m', 'function')), 'TDT not found in path, please add')
    % SPM/AFNI
    % addpath('$ADD FULL PATH TO SPM/AFNI (if you need them) AS STRING OR MAKE THIS LINE A COMMENT IF IT IS ALREADY$')
    % assert((~isempty(which('spm.m', 'function')) || ~isempty(which('BrikInfo.m', 'function'))) , 'Neither SPM nor AFNI found in path, please add (or remove this assert if you really dont need to read brain images)')


    % SET DEFAULTS
    % ------------
    cfg = decoding_defaults;


    % SET ANALYSIS
    % ------------
    % Set the analysis that should be performed (default is 'searchlight')
    cfg.analysis = 'wholebrain'; % standard alternatives: 'wholebrain', 'ROI' (pass ROIs in cfg.files.mask, see below)


    % SET MASK
    % --------
    % Set the filename of your brain mask (or your ROI masks as cell matrix) 
    % for searchlight or wholebrain e.g. 'c:\exp\glm\model_button\mask.img' OR 
    % for ROI e.g. {'c:\exp\roi\roimaskleft.img', 'c:\exp\roi\roimaskright.img'}
    % You can also use a mask file with multiple masks inside that are
    % separated by different integer values (a "multi-mask")
    % the mask is resampled ONCE, to the first subject's beta space, and reused.
    % That is only valid while every subject is in the same space, so each
    % subsequent subject is checked against the reference rather than assumed -
    % a mismatch used to pass silently and decode through a misaligned mask.
    beta1 = fullfile(firstsubjdirs{sub},'beta_0001.nii');

    if sub == 1
        mask_img = fmri_mask_image(which(mask_name),'noverbose');
        target = fmri_data(beta1,'noverbose');
        mask2write = resample_space(mask_img,target);
        mask = fullfile(firstmodelTDTdir_mask,mask_name);
        write(mask2write,'fname',mask,'overwrite');
        ref_vol = spm_vol(beta1);
        ref_vol = ref_vol(1);
    else
        this_vol = spm_vol(beta1);
        this_vol = this_vol(1);
            if ~isequal(this_vol.dim,ref_vol.dim) || ~isequal(this_vol.mat,ref_vol.mat)
                error('\n%s is not in the same space as %s, so the mask resampled to the first subject does not apply to it - please check before proceeding',firstsubjs{sub},firstsubjs{1});
            end
    end

    cfg.files.mask = mask;


    % SET OUTPUT DIR
    % --------------
    % Set the output directory where data will be saved, e.g. 'c:\exp\results\buttonpress'
    cfg.results.dir = fullfile(firstsubjdirs{sub},'TDT',analysis_name,[cfg.analysis '_' mask_name(1:end-4) '_' results_suffix]);

    if ~isfolder(cfg.results.dir)
        mkdir(cfg.results.dir);
    end


    % SET SPM.MAT PATH
    % ----------------
    % Set the filepath where your SPM.mat and all related betas are, e.g. 'c:\exp\glm\model_button'
    beta_loc = firstsubjdirs{sub};


    % SET LABEL NAMES
    % ---------------
    % Set the label names to the regressor names which you want to use for 
    % decoding, e.g. 'button left' and 'button right'
    % don't remember the names? -> run display_regressor_names(beta_loc)
    % infos on '*' (wildcard) or regexp -> help decoding_describe_data
    % indexed by POSITION in conds2include, not by the condition number itself:
    % labelnames{label} left gaps whenever conds2include was not contiguous and
    % starting at 1 (e.g. [2 4] populated {2} and {4} and left {1} and {3}
    % empty), and those empty cells went straight into decoding_describe_data
    % and every figure label
    labelnames = cell(1,numel(conds2include));
        for label = 1:numel(conds2include)
            labelnames{label} = char(spmdotmats{1,sub}.SPM.Sess(1).U(conds2include(label)).name); % labels in SPM.Sess(1).U.name, which weirdly have a blank space added to names in DSGN.conditions
        end


    % SET ADDITIONAL PARAMETERS
    % -------------------------
    % Set additional parameters manually if you want (see decoding.m or
    % decoding_defaults.m). Below some example parameters that you might want 
    % to use a searchlight with radius 12 mm that is spherical:

    % cfg.searchlight.unit = 'mm';
    % cfg.searchlight.radius = 12; % if you use this, delete the other searchlight radius row at the top!
    % cfg.searchlight.spherical = 1;
    % cfg.verbose = 2; % you want all information to be printed on screen
    % cfg.decoding.train.classification.model_parameters = '-s 0 -t 0 -c 1 -b 0 -q'; 


    % ENABLE SCALING MIN0MAX1
    % -----------------------
    % (otherwise libsvm can get VERY slow)
    % if you dont need model parameters, and if you use libsvm, use:
    % scaling per the scaling_regime option set at the top, rather than the
    % hardcoded 'all' this used to use. 'all' is TDT's own recommendation and is
    % label-blind, so it cannot carry condition information across the
    % train/test split; 'strict' estimates on the training fold only, for when
    % you would rather not estimate any quantity on data the classifier is about
    % to be tested on. Same choice, same reasoning, as in
    % LaBGAScore_decoding_SVM_between_subjects.m
    cfg.scale.method = 'min0max1';
        switch scaling_regime
            case 'kernel'
                cfg.scale.estimation = 'all';
            case 'strict'
                cfg.scale.estimation = 'across';
            otherwise
                error('\nscaling_regime must be ''kernel'' or ''strict'', not ''%s''',scaling_regime);
        end

        if sub == 1
            fprintf('\nscaling: %s/%s (regime ''%s'')\n',cfg.scale.method,cfg.scale.estimation,scaling_regime);
        end

    % if you like to change the decoding software (default: libsvm):
    % cfg.decoding.software = 'liblinear'; % for more, see decoding_toolbox\decoding_software\. 
    % Note: cfg.decoding.software and cfg.software are easy to confuse.
    % cfg.decoding.software contains the decoding software (standard: libsvm)
    % cfg.software contains the data reading software (standard: SPM/AFNI)

    % Some other cool stuff
    % Check out 
    %   combine_designs(cfg, cfg2)
    % if you like to combine multiple designs in one cfg.


    % DECIDE WHETHER YOU WANT TO SEE THE SEARCHLIGHT/ROI/... DURING DECODING
    % ----------------------------------------------------------------------
    cfg.plot_selected_voxels = 0; % 0: no plotting, 1: every step, 2: every second step, 100: every hundredth step... 0 as in LaBGAScore_decoding_SVM_between_subjects.m: under headless publishing the plots are pure cost


    % ADD ADDTIONAL OUTPUT MEASURES IF YOU LIKE
    % -----------------------------------------
    % See help decoding_transform_results for possible measures

    cfg.results.output = {'confusion_matrix'}; % 'accuracy_minus_chance' by default

    % You can also use all methods that start with "transres_", e.g. use
    %   cfg.results.output = {'SVM_pattern'};
    % will use the function transres_SVM_pattern.m to get the pattern from 
    % linear svm weights (see Haufe et al, 2015, Neuroimage)


    % NOTHING NEEDS TO BE CHANGED BELOW FOR A STANDARD LEAVE ONE-RUN OUT
    % CROSS-VALIDATED ANALYSIS
    % ------------------------------------------------------------------

    % The following function extracts all beta names and corresponding run
    % numbers from the SPM.mat
    regressor_names = design_from_spm(beta_loc);

    % Extract all information for the cfg.files structure (labels will be [1 -1] if not changed above)
    cfg = decoding_describe_data(cfg,labelnames,conds2include ,regressor_names,beta_loc);

    % This creates the leave-one-run-out cross validation design:
    cfg.design = make_design_cv(cfg); 

    % Run decoding
    results = decoding(cfg);
    results_combined{sub} = results;
    confusion_matrices = cat(3, confusion_matrices, results.confusion_matrix.output{1});


    % MAKE FIGURE
    % -----------
    figure;
    heatmap(categorical(labelnames),categorical(labelnames),results.confusion_matrix.output{1}, 'Colormap', jet);
    figtitle = firstsubjs{sub};
    set(gca,'Title',figtitle);
    plugin_set_figure_size;
    drawnow, snapnow


    % EXTRACT ACCURACY FOR EACH CONDITION FROM CONFUSION MATRIX
    % ---------------------------------------------------------
    for cond = 1:size(conds2include,2)
        group_results{cond} = [group_results{cond};results.confusion_matrix.output{1}(cond,cond)]; 
    end

end % for loop over subjects


%% SAVE RESULTS
% -------------------------------------------------------------------------

fprintf('\n\n');
printhdr('SAVING CROSS-CLASSIFICATION ACCURACY RESULTS');
fprintf('\n\n');

savefilenamedata = fullfile(secondmodelTDTanalysisdir, [cfg.analysis '_' mask_name(1:end-4) '_' results_suffix '.mat']);
save(savefilenamedata, 'confusion_matrices', 'group_results','labelnames','-v7.3');

fprintf('\nSaved results for %s\n', analysis_name);
fprintf('\nFilename: %s\n', savefilenamedata);


%% CALCULATE AVERAGE AND DO GROUP LEVEL STATS
% -------------------------------------------------------------------------

mean_matrix = mean(confusion_matrices,3);
std_matrix = std(confusion_matrices,[],3);
median_matrix = median(confusion_matrices,3);
iqr_matrix = iqr(confusion_matrices,3);

figure;
heatmap(categorical(labelnames),categorical(labelnames),median_matrix, 'Colormap', jet);
figtitle = 'median classification accuracy';
set(gca,'Title',figtitle);
plugin_set_figure_size; % as elsewhere in the repo, rather than WindowState maximized, which sizes from the screen rather than for publishing
drawnow, snapnow

figure;
heatmap(categorical(labelnames),categorical(labelnames),iqr_matrix, 'Colormap', jet);
figtitle = 'iqr classification accuracy';
set(gca,'Title',figtitle);
plugin_set_figure_size;
drawnow, snapnow

group_results = cell2mat(group_results);
group_results_tbl = array2table(group_results);
varnames = group_results_tbl.Properties.VariableNames;

    for var = 1:size(varnames,2)
        group_results_tbl.Properties.VariableNames{var} = ['c' num2str(var)];
        group_results_tbl.Properties.VariableDescriptions{var} = labelnames{var};
    end

conds = table((1:numel(labelnames))','VariableNames',{'Conditions'});

rm = fitrm(group_results_tbl,[group_results_tbl.Properties.VariableNames{1} '-' group_results_tbl.Properties.VariableNames{end} ' ~ 1'],'WithinDesign',conds);
ranovatbl = ranova(rm);
margmeanstbl = margmean(rm,'Conditions');
posthoctbl = multcompare(rm,'Conditions'); % default Tukey-Kramer corrected

% chance level depends on the SCALE of TDT's confusion_matrix output. It is
% written in percent, so chance is 100/nConditions - but that assumption was
% silent here, and a proportion-scaled matrix would have been tested against a
% chance level 100x too high, returning a "significant" result for every
% condition. Derive the scale from the data, and say which one was used.
n_conds = numel(conds2include);

    if max(confusion_matrices(:)) <= 1.5
        chance_level = 1/n_conds;
        warning('\nconfusion matrix values look like PROPORTIONS (max %.3f), testing against chance = %.4f - please verify this is what TDT returned',max(confusion_matrices(:)),chance_level);
    else
        chance_level = 100/n_conds;
        fprintf('\nconfusion matrix in percent (max %.1f), testing against chance = %.2f%%\n',max(confusion_matrices(:)),chance_level);
    end

inputs = cell(1,numel(varnames));
p      = cell(1,numel(varnames));
stats  = cell(1,numel(varnames));

    for var = 1:size(varnames,2)
        input = confusion_matrices(var,var,:);
        inputs{var} = input(:);
        [p{var},~,stats{var}] = signrank(inputs{var},chance_level); % the missing closing parenthesis here made the whole file unparseable
        clear input;
    end
    
inputs_matrix = cell2mat(inputs);

figure;
boxplot(inputs_matrix,'Notch','on','Labels',labelnames,'Colors',lines(n_conds)); % lines(n) rather than a hardcoded 'rgbm', which only ever coloured four conditions
xlabel('condition');
    if chance_level > 1
        ylabel('% accuracy');
    else
        ylabel('proportion correct');
    end
yline(chance_level,'--k','chance');
plugin_set_figure_size;
drawnow, snapnow


%% SAVE SECOND LEVEL RESULTS
% -------------------------------------------------------------------------

fprintf('\n\n');
printhdr('SAVING CROSS-CLASSIFICATION ACCURACY GROUP LEVEL STATS');
fprintf('\n\n');

savefilenamedata = fullfile(secondmodelTDTanalysisdir, ['group_level_stats_' cfg.analysis '_' mask_name(1:end-4) '_' results_suffix '.mat']);
save(savefilenamedata, 'mean_matrix', 'median_matrix', 'std_matrix', 'iqr_matrix', 'group_results_tbl', 'rm', 'inputs_matrix', ...
    'ranovatbl', 'margmeanstbl', 'posthoctbl', 'p', 'stats', 'chance_level', 'labelnames', '-v7.3'); % the ranova, post-hoc and sign-rank results were computed but never saved

fprintf('\nSaved results for %s\n', analysis_name);
fprintf('\nFilename: %s\n', savefilenamedata);