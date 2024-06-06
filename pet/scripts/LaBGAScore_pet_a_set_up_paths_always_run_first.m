%% LaBGAScore_pet_a_set_up_paths_always_run_first.m
%
% 
% USAGE
% 
% Always run this first before you run other LaBGAScore pet second level batch scripts.
%
%
% CANLAB NOTES
%
% - standard folders and variable names are created by these scripts
%
% - in "prep_" scripts: 
%   image names, conditions, contrasts, colors, global gray/white/CSF
%   values are saved automatically in a DAT structure
% 
% - extracted fmri_data objects are saved in DATA_OBJ variables
% - contrasts are estimated and saved in DATA_OBJ_CON variables
%
% - files with these variables are saved and loaded automatically when you
%   run the scripts
%   meta-data saved in image_names_and_setup.mat
%   image data saved in data_objects.mat
%
% - you only need to run the prep_ scripts once.  After that, use 
%   b_reload_saved_matfiles.m to re-load saved files
% 
% - when all scripts working properly, run z_batch_publish_analyses.m
%   to create html report.  customize by editing z_batch_list_to_publish.m
%
% - saved in results folder:
%   figures
%   html report with figures and stats, in "published_output"
%
%
% LaBGAS NOTES
%
% - script to be run from rootdir of superdataset for your study
% - DO NOT FORGET TO MAKE STUDY-SPECIFIC CHANGES INDICATED BELOW
%
%__________________________________________________________________________
%
% modified for pet data by: Lukas Van Oudenhove
% date: KU Leuven, June 2024    
%
%__________________________________________________________________________
% @(#)% LaBGAScore_pet_a_set_up_paths_always_run_first.m         v1.0
% last modified: 2024/06/05


%% RUN PREP AND FIRST LEVEL DESIGN SCRIPT
% -------------------------------------------------------------------------

% check whether LaBGAScore_prep_s0_define_directories has been run
% STUDY-SPECIFIC: replace LaBGAScore with study name in code below

if ~exist('rootdir','var') || ~exist('githubrootdir','var')
    fprintf('\n')
    warning('rootdir and/or githubrootdir variable not found in Matlab workspace, running LaBGAScore_prep_s0_define_directories before proceeding')
    fprintf('\n')
    LaBGAScore_prep_s0_define_directories;
    cd(rootdir);
else
    cd(rootdir);
end

% specify tracer name

modelname = 'trc-DPA714'; % give your model the name of your tracer used in LaBGAScore_pet_model script


%% SET DEFAULT USER OPTIONS
% -------------------------------------------------------------------------

% STUDY-SPECIFIC: add study name and model name to script name

LaBGAScore_pet_a2_set_default_options;

    
%% MAKE SURE DEPENDENCIES ARE ON MATLAB PATH
% -------------------------------------------------------------------------

% check whether spm subdirs are on path, add if needed

spmcanonicaldir = fullfile(spmrootdir,'canonical');
    if sum(contains(matlabpath,spmcanonicaldir)) == 0
        addpath(spmcanonicaldir,'-end');
        warning('\nadding %s to end of Matlab path',spmcanonicaldir)
    end
spmconfigdir = fullfile(spmrootdir,'config');
    if sum(contains(matlabpath,spmconfigdir)) == 0
        addpath(spmconfigdir,'-end');
        warning('\nadding %s to end of Matlab path',spmconfigdir)
    end
spmmatlabbatchdir = fullfile(spmrootdir,'matlabbatch');
    if sum(contains(matlabpath,spmmatlabbatchdir)) == 0
        addpath(spmmatlabbatchdir,'-end');
        warning('\nadding %s to end of Matlab path',spmmatlabbatchdir)
    end
spmtoolboxdir = fullfile(spmrootdir,'toolbox');
    if sum(contains(matlabpath,spmtoolboxdir)) == 0
        addpath(spmtoolboxdir,'-end');
        warning('\nadding %s to end of Matlab path',spmtoolboxdir)
    end
    
% check whether CANlab Github repos are cloned and on Matlab path, clone and/or add if needed

  % CANLABCORE
    canlabcoredir = fullfile(githubrootdir,'CanlabCore');
        if ~isfolder(canlabcoredir) % canlabcore not yet cloned
          canlabcoreurl = "https://github.com/canlab/CanlabCore.git";
          canlabcoreclonecmd = ['git clone ' canlabcoreurl];
          cd(githubrootdir);
          [status,cmdout] = system(canlabcoreclonecmd);
          disp(cmdout);
              if status == -0
                  addpath(genpath(canlabcoredir,'-end'));
                  warning('\ngit succesfully cloned %s to %s and added repo to Matlab path',canlabcoreurl, canlabcoredir)
              else
                  error('\ncloning %s into %s failed, please try %s in linux terminal before proceeding, or use Gitkraken',canlabcoreurl,canlabcoredir,canlabcoreclonecmd)
              end
          cd(rootdir);
          clear status cmdout
        elseif ~exist('fmri_data.m','file') % canlabcore cloned but not yet on Matlab path
            addpath(genpath(canlabcoredir),'-end');
        end
        
  % CANLABPRIVATE
    canlabprivdir = fullfile(githubrootdir,'CanlabPrivate');
        if ~isfolder(canlabprivdir) % canlabprivate not yet cloned
          canlabprivurl = "https://github.com/canlab/CanlabPrivate.git";
          canlabprivclonecmd = ['git clone ' canlabprivurl];
          cd(githubrootdir);
          [status,cmdout] = system(canlabprivclonecmd);
          disp(cmdout);
              if status == -0
                  addpath(genpath(canlabprivdir,'-end'));
                  warning('\ngit succesfully cloned %s to %s and added repo to Matlab path',canlabprivurl, canlabprivdir)
              else
                  error('\ncloning %s into %s failed, please try %s in linux terminal before proceeding, or use Gitkraken',canlabprivurl,canlabprivdir,canlabprivclonecmd)
              end
          cd(rootdir);
          clear status cmdout
        elseif ~exist('power_calc.m','file') % canlabprivate cloned but not yet on Matlab path
            addpath(genpath(canlabprivdir),'-end');
        end
        
  % CANLAB HELP EXAMPLES (LaBGAS fork)
    canlabhelpdir = fullfile(githubrootdir,'CANlab_help_examples');
        if ~isfolder(canlabhelpdir) % CANlab_help_examples not yet cloned
          canlabhelpurl = "https://github.com/labgas/CANlab_help_examples.git";
          canlabhelpclonecmd = ['git clone ' canlabhelpurl];
          cd(githubrootdir);
          [status,cmdout] = system(canlabhelpclonecmd);
          disp(cmdout);
              if status == -0
                  addpath(genpath(canlabhelpdir,'-end'));
                  warning('\ngit succesfully cloned %s to %s and added repo to Matlab path',canlabhelpurl, canlabhelpdir)
              else
                  error('\ncloning %s into %s failed, please try %s in linux terminal before proceeding, or use Gitkraken',canlabhelpurl,canlabhelpdir,canlabhelpclonecmd)
              end
          cd(rootdir);
          clear status cmdout
        elseif ~exist('a0_begin_here_readme.m','file') % CANlab_help_examples cloned but not yet on Matlab path
            addpath(genpath(canlabhelpdir),'-end');
        end
        
  % CANLAB SINGLE TRIALS
        canlabsingletrialsdir = fullfile(githubrootdir,'canlab_single_trials');
        if ~isfolder(canlabsingletrialsdir) % canlab_single_trials not yet cloned
          canlabsingletrialsurl = "https://github.com/labgas/canlab_single_trials.git";
          canlabsingletrialsclonecmd = ['git clone ' canlabsingletrialsurl];
          cd(githubrootdir);
          [status,cmdout] = system(canlabsingletrialsclonecmd);
          disp(cmdout);
              if status == -0
                  addpath(genpath(canlabsingletrialsdir,'-end'));
                  warning('\ngit succesfully cloned %s to %s and added repo to Matlab path',canlabsingletrialsurl, canlabsingletrialsdir)
              else
                  error('\ncloning %s into %s failed, please try %s in linux terminal before proceeding, or use Gitkraken',canlabsingletrialsurl,canlabsingletrialsdir,canlabsingletrialsclonecmd)
              end
          cd(rootdir);
          clear status cmdout
        elseif ~exist('fmri_data_st.m','file') % canlab_single_trials cloned but not yet on Matlab path
            addpath(genpath(canlabsingletrialsdir),'-end');
        end

        
%% SET BASE DIRECTORY AND CREATE STANDARD SUBDIR STRUCTURE
% -------------------------------------------------------------------------

% Pet directory in derivatives subdataset

derivrootdir = fileparts(derivdir);
derivpetdir = fullfile(derivrootdir,['pet_' modelname]);

% Pet directory in firstlevel subdataset

datadir = fullfile(rootdir,'firstlevel',['pet_' modelname]);
    if ~exist(datadir, 'dir')
        error('\nfirstleveldir for tracer %s does not exist, please check naming and consistency with %s, or run LaBGAScore_pet_model script first', modelname, fullfile(rootdir,'firstlevel')) 
    end

% Base directory for second level model with standard subdirs

basedir = fullfile(rootdir,'secondlevel',['pet_' modelname]);
    if ~exist(basedir, 'dir')
        mkdir(basedir); 
    end
maskdir = fullfile(basedir,'masks');
    if ~exist(maskdir, 'dir')
        mkdir(maskdir); 
    end
roidir = fullfile(maskdir,'rois');
    if ~exist(maskdir, 'dir')
        mkdir(maskdir); 
    end
addpath(genpath(maskdir),'-end');
resultsdir = fullfile(basedir, 'results');
    if ~exist(resultsdir, 'dir')
        mkdir(resultsdir); 
    end
figsavedir = fullfile(resultsdir, 'figures');
    if ~exist(figsavedir, 'dir')
        mkdir(figsavedir); 
    end
notesdir = fullfile(resultsdir, 'notes');
    if ~exist(notesdir, 'dir')
        mkdir(notesdir); 
    end
htmlsavedir = fullfile(resultsdir,'html');
    if ~exist(htmlsavedir,'dir')
        mkdir(htmlsavedir);
    end

    
%% DEFINE HELPER FUNCTION CALLED BY LATER SCRIPTS
% -------------------------------------------------------------------------

dashes = '----------------------------------------------';
printstr = @(dashes) disp(dashes);
printhdr = @(str) fprintf('%s\n%s\n%s\n', dashes, str, dashes);