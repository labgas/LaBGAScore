%% LaBGAScore_prep_s0_define_directories
%
% This script will define the paths for the standard BIDS-compliant
% directory structure for (datalad) datasets for LaBGAS neuroimaging projects
% 
% USAGE
% Script should be run from the root directory of your superdataset, e.g.
% /data/proj_discoverie, on Linux OS
%
% DEPENDENCIES
% spm12 added to your Matlab path WITHOUT subdirectories
% no spm functions are called by this script, 
% but spmrootdir is defined automatically, for use in later scripts
%
%__________________________________________________________________________
%
% author: Lukas Van Oudenhove
% date:   November, 2021
%
%__________________________________________________________________________
% @(#)% LaBGAScore_prep_s0_define_directories.m         v1.1       
% last modified: 2022/03/14
% 
%
%% DEFINE DIRECTORIES AND ADD CODE DIR TO MATLAB PATH
%--------------------------------------------------------------------------
rootdir = pwd;
sourcedir = fullfile(rootdir,'sourcedata');
BIDSdir = fullfile(rootdir,'BIDS');
codedir = fullfile(rootdir,'code');
derivdir = fullfile(rootdir,'derivatives','fmriprep');
matlabpath = path;

    if ~exist('spm.m','file')
        spmpathcommand = "addpath('your_spm_rootdir','-end')";
        error('spm12 not found on Matlab path, please add WITHOUT subfolders using the Matlab GUI or type %s in Matlab terminal before proceeding',spmpathcommand)
    else
        spmrootdir = which('spm.m');
        spmrootdir = strsplit(spmrootdir,'/spm.m');
        spmrootdir = spmrootdir{1,1};
    end

if sum(contains(matlabpath,codedir)) == 0
    addpath(genpath(codedir),'-end');
    warning('adding %s to end of Matlab path',codedir)
end

%% READ IN SUBJECT LISTS AND COMPARE THEM
%--------------------------------------------------------------------------
sourcelist = dir(fullfile(sourcedir,'sub-*'));
sourcesubjs = cellstr(char(sourcelist(:).name));
BIDSlist = dir(fullfile(BIDSdir,'sub-*'));
BIDSsubjs = cellstr(char(BIDSlist(:).name));
derivlist = dir(fullfile(derivdir,'sub-*'));
derivlist = derivlist([derivlist(:).isdir]);
derivsubjs = cellstr(char(derivlist.name));

if isequal(sourcesubjs,BIDSsubjs,derivsubjs)
    warning('\n numbers and names of subjects in %s, %s, and %s match - good to go',sourcedir,BIDSdir,derivdir);
else
    error('\n numbers and names of subjects in %s, %s, and %s do not match - please check before proceeding and make sure your file organization is consistent with LaBGAS conventions',sourcedir,BIDSdir,derivdir);
end


%% CREATE CELL ARRAYS WITH FULL PATHS FOR SUBJECT DIRECTORIES
%--------------------------------------------------------------------------
for sourcesub = 1:size(sourcesubjs,1)
    sourcesubjdirs{sourcesub,1} = fullfile(sourcelist(sourcesub).folder,sourcelist(sourcesub).name);
end

for BIDSsub = 1:size(BIDSsubjs,1)
    BIDSsubjdirs{BIDSsub,1} = fullfile(BIDSlist(BIDSsub).folder,BIDSlist(BIDSsub).name);
end

for derivsub = 1:size(derivsubjs,1)
    derivsubjdirs{derivsub,1} = fullfile(derivlist(derivsub).folder,derivlist(derivsub).name);
end


%% CLEAN UP OBSOLETE VARIABLES
%--------------------------------------------------------------------------
clear sourcelist BIDSlist derivlist sourcesub BIDSsub derivsub 