%% LaBGAScore_prep_s2_smooth.m
%
%
% *USAGE*
%
% This script smooths fMRIprep output images for a single-session dataset,
% more specifically
%
% 1. define directories and subject lists by calling
%   <study_prefix>_prep_s0_define_directories
%
% 2. check the subjects requested in subjs2smooth against the subjects
%   present in derivdir, and error out on any that are missing
%
% 3. for each subject, in derivatives/fmriprep/<sub>/func
%   a) gunzip the *preproc_bold*.nii.gz images
%   b) build and save an spm12 smoothing batch as <sub>_smooth.mat
%   c) run the batch, producing <prefix>*.nii
%   d) gzip the smoothed images
%   e) delete all unzipped .nii images again
%   NOTE: spm functions called by this script
%       a) spm_select - to collect the unzipped images
%       b) spm_jobman - to run the smoothing batch
%   NOTE: the batch is saved before it is run, so the smoothing is
%       reproducible from the batch file alone
%
% Script should be run from the root directory of the superdataset, e.g.
% /data/proj_discoverie
% The script is generic, i.e. it does not require study-specific adaptions,
% but you can change some default options if required.
% For datasets with more than one session, use
% LaBGAScore_prep_s2_smooth_multisess.m instead.
%
%
% *OPTIONS*
%
% * study_prefix    prefix used for all scripts of a given study, STUDY-SPECIFIC
%
% * fwhm            smoothing kernel width in mm
%
% * prefix          string defining prefix of choice for smoothing images
%
% * subjs2smooth    cell array of subjects in derivdir you want to smooth, empty cell array (default) if you want to loop over all subjects
%
%
% *DEPENDENCIES*
%
% 1. LaBGAScore Github repo on Matlab path, with subfolders
%   https://github.com/labgas/LaBGAScore
% 2. spm12 on Matlab path, without subfolders
%   will be checked by calling LaBGAScore_prep_s0_define_directories
%
%
% *NOTES*
%
% INPUTS: preprocessed .nii.gz images outputted by fMRIprep; variables
% created by running LaBGAScore_prep_s0_define_directories from the root
% directory of your (super)dataset
%
% OUTPUT: smoothed .nii.gz images, plus the spm batch that produced them
%
% prefix and fwhm are independent: the prefix is a label, not something
% derived from the kernel, so keeping them consistent is up to you. The
% firstlevel DSGN.funcnames must glob for the same prefix (every shipped
% example uses 's6*').
%
% -------------------------------------------------------------------------
%
% modified by: Lukas Van Oudenhove
%
% date:   November, 2021
%
% -------------------------------------------------------------------------
%
% LaBGAScore_prep_s2_smooth.m         v1.5
%
% last modified: 2026/09/30
%
%
%% SET SMOOTHING OPTIONS, AND SUBJECTS
%--------------------------------------------------------------------------

study_prefix = '';  % STUDY-SPECIFIC
fwhm = 6;           % kernel width in mm
prefix = 's6-';     % prefix for name of smoothed images
subjs2smooth = {};  % enter subjects separated by comma if you only want to smooth selected subjects e.g. {'sub-01','sub-02'}


%% DEFINE DIRECTORIES
%--------------------------------------------------------------------------

eval([study_prefix '_prep_s0_define_directories']);


%% UNZIP IMAGES, SMOOTH, ZIP, SMOOTHED IMAGES, AND DELETE ALL UNZIPPED IMAGES
%----------------------------------------------------------------------------

% resolve WHICH subjects to loop over first, so the body below exists in a
% single copy - it used to be duplicated verbatim between the subjs2smooth
% and the all-subjects branch, which meant every fix had to be applied twice
if ~isempty(subjs2smooth)
    [C,ia,~] = intersect(derivsubjs,subjs2smooth);
        if ~isequal(C,subjs2smooth')
            error('\n subject %s defined in subjs2smooth not present in %s, please check before proceeding',subjs2smooth{~ismember(subjs2smooth,C)},derivdir);
        end
    subjidx = ia';
else
    subjidx = 1:size(derivsubjdirs,1);
end

for sub = subjidx

    cd(fullfile(derivsubjdirs{sub,:},'func'));

    % unzip .nii.gz files
    % gunzip RETURNS the files it created; using that list rather than a
    % wildcard keeps both the smoothing below and the cleanup afterwards
    % scoped to the images this iteration actually unzipped, instead of every
    % .nii that happens to sit in this directory
    unzipped = gunzip('*preproc_bold*.nii.gz');

    % write smoothing spm batch
    % spm_select is called per unzipped file so that 4D volumes are still
    % frame-expanded (ExtFPList), while the selection stays restricted to the
    % files just unzipped
    clear matlabbatch;
    scans = {};
        for u = 1:numel(unzipped)
            [~,uname,uext] = fileparts(unzipped{u});
            scans = [scans; cellstr(spm_select('ExtFPList',pwd,['^' regexptranslate('escape',[uname uext]) '$'],Inf))]; %#ok<AGROW>
        end
    kernel = ones(1,3).*fwhm;
    matlabbatch{1}.spm.spatial.smooth.data = scans;
    matlabbatch{1}.spm.spatial.smooth.fwhm = kernel;
    matlabbatch{1}.spm.spatial.smooth.dtype = 0;
    matlabbatch{1}.spm.spatial.smooth.im = 0;
    matlabbatch{1}.spm.spatial.smooth.prefix = prefix;

    % save batch and run
    save(fullfile(pwd,[derivsubjs{sub,:} '_smooth.mat']),'matlabbatch');
    spm_jobman('initcfg');
    spm_jobman('run',matlabbatch);

    % zip the smoothed images, then delete every image this iteration
    % unzipped or created, addressing both BY NAME rather than by wildcard.
    % The old gzip('s6*') hardcoded the DEFAULT prefix while prefix is an
    % option, so changing prefix matched nothing, left the smoothed images
    % unzipped, and the following delete('*.nii') then removed them - i.e. it
    % silently destroyed the smoothing output
        for u = 1:numel(unzipped)
            [udir,uname,uext] = fileparts(unzipped{u});
            smoothedfile = fullfile(udir,[prefix uname uext]);
                if isfile(smoothedfile)
                    gzip(smoothedfile);
                    delete(smoothedfile);
                else
                    warning('\nexpected smoothed image %s not found, so it was neither zipped nor deleted - please check',smoothedfile);
                end
            delete(unzipped{u});
        end

end % for loop over subjects

% return to the superdataset root: s0 sets rootdir = pwd, so leaving pwd on
% the last subject would make any script run afterwards in the same Matlab
% session derive its paths from the wrong place
cd(rootdir);
