%% LaBGAScore_prep_s0_define directories
%
% This script will define the paths for the standard BIDS-compliant
% directory structure for (datalad) datasets for LaBGAS neuroimaging projects
% 
% USAGE
% Script should be run from the root directory of your superdataset, e.g.
% /data/proj_discoverie, on Linux OS
%
%__________________________________________________________________________
%
% author: Lukas Van Oudenhove
% date:   November, 2021
%
%__________________________________________________________________________
% 
%
%% DEFINE DIRECTORIES AND ADD CODE DIR TO MATLAB PATH
%--------------------------------------------------------------------------
rootdir = pwd;
sourcedir = fullfile(rootdir,'sourcedata');
BIDSdir = fullfile(rootdir,'BIDS');
codedir = fullfile(rootdir,'code');
derivdir = fullfile(rootdir,'derivatives','fmriprep');

addpath(genpath(codedir));


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
    disp('numbers or names of subjects in sourcedata, BIDS, and derivatives/fmriprep directory match - good to go');
else
    warning('numbers or names of subjects in sourcedata, BIDS, and derivatives/fmriprep directory do not match - please check');;
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