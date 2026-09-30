%% LaBGAScore_prep_s0_define_directories.m
%
%
% *USAGE*
%
% This script defines the standard LaBGAS BIDS-compliant directory structure
% and subject lists for a neuroimaging (datalad) dataset, more specifically
%
% 1. set the study prefix that binds all of a study's scripts together
%   NOTE: study_prefix is STUDY-SPECIFIC and has no usable default - the
%       other prep scripts reach this one through
%       eval([study_prefix '_prep_s0_define_directories']), which is not a
%       valid Matlab identifier when study_prefix is left empty
%
% 2. define the directories of the standard LaBGAS layout: rootdir (= pwd),
%   sourcedir, BIDSdir, codedir, derivdir (derivatives/fmriprep), and
%   githubrootdir
%
% 3. check that spm12 is on the Matlab path, and derive spmrootdir from it
%   NOTE: no spm function is called by this script; spmrootdir is defined
%       for use by later scripts, particularly firstlevel s2
%   NOTE: errors out with the exact addpath command to run if spm.m is not
%       found on the path
%
% 4. add the study's code dir to the Matlab path if not already present
%
% 5. read the sub-* lists from sourcedir, BIDSdir and derivdir and compare them
%   NOTE: a mismatch is an ERROR rather than a warning, because every later
%       script indexes sourcesubjs/BIDSsubjs/derivsubjs and their *subjdirs
%       counterparts POSITIONALLY - sourcesubjdirs{k}, BIDSsubjdirs{k} and
%       derivsubjdirs{k} must be the same participant
%
% 6. build cell arrays of full paths for the subject dirs of all three trees
%
% This script should be run from the root directory of the superdataset, e.g.
% /data/proj_discoverie, since rootdir is set to pwd - starting anywhere else
% silently points the whole pipeline at the wrong tree.
% It can be used standalone, but is typically called from the subsequent
% scripts in the standard LaBGAS workflow, from firstlevel, and from the
% second level templates in the CANlab_help_examples LaBGAS fork.
% The directory-detection logic itself is generic and does not require
% study-specific adaptions, beyond setting study_prefix below.
%
%
% *OPTIONS*
%
% * study_prefix    prefix used for all scripts of a given study, STUDY-SPECIFIC
%
%
% *DEPENDENCIES*
%
% spm12 on Matlab path WITHOUT subdirectories
% no spm functions are called by this script,
% but spmrootdir is defined automatically, for use in later scripts
%
%
% *NOTES*
%
% OUTPUTS (in the base workspace): rootdir, githubrootdir, sourcedir, BIDSdir,
% codedir, derivdir, spmrootdir, sourcesubjs/BIDSsubjs/derivsubjs, and
% sourcesubjdirs/BIDSsubjdirs/derivsubjdirs
%
% githubrootdir is hardcoded to '/data/master_github_repos', i.e. it
% assumes all LaBGAS Github repos are cloned locally under that path;
% adapt this line if your local repo location differs
%
% See prep/README.md for how this script hands off to the rest of prep and to
% firstlevel
%
% -------------------------------------------------------------------------
%
% modified by: Lukas Van Oudenhove
%
% date:   November, 2021
%
% -------------------------------------------------------------------------
%
% LaBGAScore_prep_s0_define_directories.m         v1.4
%
% last modified: 2026/09/30
%
%
%% SET STUDY PREFIX FOR USE IN ALL SUBSEQUENT SCRIPTS
%--------------------------------------------------------------------------
study_prefix = ''; % STUDY-SPECIFIC


%% DEFINE DIRECTORIES AND ADD CODE DIR TO MATLAB PATH
%--------------------------------------------------------------------------
rootdir = pwd;
githubrootdir = '/data/master_github_repos'; %dir where all your Github repos are cloned locally
sourcedir = fullfile(rootdir,'sourcedata');
BIDSdir = fullfile(rootdir,'BIDS');
codedir = fullfile(rootdir,'code');
derivdir = fullfile(rootdir,'derivatives','fmriprep');
matlabpath = path;

    if ~exist('spm.m','file')
        spmpathcommand = "addpath('your_spm_rootdir','-end')";
        error('\nspm12 not found on Matlab path, please add WITHOUT subfolders using the Matlab GUI or type %s in Matlab terminal before proceeding',spmpathcommand)
    else
        spmrootdir = fileparts(which('spm.m')); % fileparts rather than strsplit on '/spm.m', which hardcoded a forward slash and hence did not work on Windows
    end

if sum(contains(matlabpath,codedir)) == 0
    addpath(genpath(codedir),'-end');
    warning('\nadding %s to end of Matlab path',codedir)
end

%% READ IN SUBJECT LISTS AND COMPARE THEM
%--------------------------------------------------------------------------
% all three listings are filtered to DIRECTORIES: without this a stray file
% named sub-* in sourcedata or BIDS is counted as a subject and trips the
% three-way check below
sourcelist = dir(fullfile(sourcedir,'sub-*'));
sourcelist = sourcelist([sourcelist(:).isdir]);
sourcesubjs = cellstr(char(sourcelist(:).name));
BIDSlist = dir(fullfile(BIDSdir,'sub-*'));
BIDSlist = BIDSlist([BIDSlist(:).isdir]);
BIDSsubjs = cellstr(char(BIDSlist(:).name));
derivlist = dir(fullfile(derivdir,'sub-*'));
derivlist = derivlist([derivlist(:).isdir]);
derivsubjs = cellstr(char(derivlist.name));

% an empty listing has to be reported as such: cellstr(char([])) yields {''}
% rather than an empty cell, so "no subjects found" would otherwise surface
% below as a confusing list mismatch
if isempty(sourcelist) || isempty(BIDSlist) || isempty(derivlist)
    error('\nno sub-* directories found in one or more of %s, %s, %s - please check that your dataset is organized according to LaBGAS convention, and that fMRIprep has been run',sourcedir,BIDSdir,derivdir);
end

if isequal(sourcesubjs,BIDSsubjs,derivsubjs)
    warning('\nnumbers and names of subjects in %s, %s, and %s match - good to go',sourcedir,BIDSdir,derivdir);
else
    error('\nnumbers and names of subjects in %s, %s, and %s do not match - please check before proceeding and make sure your file organization is consistent with LaBGAS conventions',sourcedir,BIDSdir,derivdir);
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
% the loop counters and the path snapshot are of no use to calling scripts,
% and leaving them in the base workspace invites a later script to reuse one
% by accident
clear sourcesub BIDSsub derivsub matlabpath spmpathcommand
clear sourcelist BIDSlist derivlist sourcesub BIDSsub derivsub 