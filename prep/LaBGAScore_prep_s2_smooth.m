%% LaBGAScore_prep_s2_smooth
%
% This script unzips fMRIprep output images, smooth them, zip the
% smoothed images, and delete all the unzipped images again
%
% USAGE
%
% Script should be run from the root directory of the superdataset, e.g.
% /data/proj_discoverie
% The script is generic, i.e. it does not require study-specific adaptions,
% but you can change some default options if required
%
%
% DEPENDENCIES
%
% 1. LaBGAScore Github repo on Matlab path, with subfolders
%   https://github.com/labgas/LaBGAScore
% 2. spm12 on Matlab path, without subfolders
%   will be checked by calling LaBGAScore_prep_s0_define_directories
%
%
% INPUTS
%
% preprocessed .nii.gz images outputted by fMRIprep
% variables created by running LaBGAScore_prep_s0_define_directories from
% the root directory of your (super)dataset
%
%
% OUTPUT
%
% smoothed .nii.gz images
%
%
% OPTIONS
%
% 1. fwhm
%   smoothing kernel width in mm
% 2. prefix
%   string defining prefix of choice for smoothing images
% 3. subjs2smooth
%   cell array of subjects in derivdir you want to smooth, empty cell array
%   if you want to loop over all subjects
%
%__________________________________________________________________________
%
% author: Lukas Van Oudenhove
% date:   November, 2021
%
%__________________________________________________________________________
% @(#)% LaBGAScore_prep_s2_smooth.m         v1.1       
% last modified: 2022/03/15


%% SET SMOOTHING OPTIONS, AND SUBJECTS
%--------------------------------------------------------------------------

fwhm = 6; % kernel width in mm
prefix = 's6-'; % prefix for name of smoothed images
subjs2smooth = {}; % enter subjects separated by comma if you only want to smooth selected subjects e.g. {'sub-01','sub-02'}


%% DEFINE DIRECTORIES
%--------------------------------------------------------------------------

LaBGAScore_prep_s0_define_directories;


%% UNZIP IMAGES, SMOOTH, ZIP, SMOOTHED IMAGES, AND DELETE ALL UNZIPPED IMAGES
%----------------------------------------------------------------------------

if ~isempty(subjs2smooth)
    [C,ia,~] = intersect(derivsubjs,subjs2smooth);
    
    if ~isequal(C,subjs2smooth)
        error('\n subject %s defined in subjs2smooth not present in %s, please check before proceeding',subjs2smooth{~ismember(subjs2smooth,C)},derivdir);
    else
        
        for sub=ia'
            cd([derivsubjdirs{sub,:},'/func']);
            % unzip .nii.gz files
            gunzip('*preproc_bold*.nii.gz');
            % write smoothing spm batch
            clear matlabbatch;
            matlabbatch = struct([]);
            scans=spm_select('ExtFPList',pwd,'.*\.nii$',Inf);
            kernel = ones(1,3).*fwhm;
            matlabbatch{1}.spm.spatial.smooth.data = cellstr(scans);
            matlabbatch{1}.spm.spatial.smooth.fwhm = kernel;
            matlabbatch{1}.spm.spatial.smooth.dtype = 0;
            matlabbatch{1}.spm.spatial.smooth.im = 0;
            matlabbatch{1}.spm.spatial.smooth.prefix = prefix;
            % save batch and run
            eval(['save ' derivsubjs{sub,:} '_smooth.mat matlabbatch']); 
            spm_jobman('initcfg');
            spm_jobman('run',matlabbatch);
            % zip smoothed files
            gzip('s6*');
            % delete all unzipped files
            delete('*.nii');
        end % for loop over subjs2smooth
        
    end % if loop checking intersection of subjs2smooth and subjdirs
    
else
    
    for sub=1:size(derivsubjdirs,1)
        cd([derivsubjdirs{sub,:},'/func']);
        % unzip .nii.gz files
        gunzip('*preproc_bold*.nii.gz');
        % write smoothing spm batch
        clear matlabbatch;
        matlabbatch = struct([]);
        scans=spm_select('ExtFPList',pwd,'.*\.nii$',Inf);
        kernel = ones(1,3).*fwhm;
        matlabbatch{1}.spm.spatial.smooth.data = cellstr(scans);
        matlabbatch{1}.spm.spatial.smooth.fwhm = kernel;
        matlabbatch{1}.spm.spatial.smooth.dtype = 0;
        matlabbatch{1}.spm.spatial.smooth.im = 0;
        matlabbatch{1}.spm.spatial.smooth.prefix = prefix;
        % save batch and run
        eval(['save ' derivsubjs{sub,:} '_smooth.mat matlabbatch']); 
        spm_jobman('initcfg');
        spm_jobman('run',matlabbatch);
        % zip smoothed files
        gzip('s6*');
        % delete all unzipped files
        delete('*.nii');
    end % for loop over subjdirs
    
end % if loop checking smoothing option

cd(rootdir);