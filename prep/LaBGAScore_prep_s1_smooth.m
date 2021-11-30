%% LaBGAScore_prep_s1_smooth
%
% This script will unzip fMRIprep output images, smooth them, zip the
% smoothed images, and delete all the unzipped images again
%
% IMPORTANT NOTE
% Run this script from the root dir of your (super)dataset as it first
% calls LaBGAScore_prep_s0_define_directories which uses a relative path
% 
% DEPENDENCIES
% SPM12 on your Matlab path
% 
% INPUTS
% preprocessed .nii.gz images outputted by fMRIprep
% variables created by running LaBGAScore_prep_s0_define_directories from
% the root directory of your (super)dataset
%
% OUTPUT
% smoothed .nii.gz images
%
% OPTIONS
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
% @(#)% LaBGAScore_prep_s1_smooth.m         v1.0       
% last modified: 2021/06/04
%
%
%% SET SMOOTHING OPTIONS, AND SUBJECTS
%--------------------------------------------------------------------------

fwhm = 6; % kernel width in mm
prefix = 's6-'; % prefix for name of smoothed images
subjs2smooth = {}; % specify if you only want to smooth a subset of all subjects in derivdir, otherwise leave cell array empty


%% DEFINE DIRECTORIES
%--------------------------------------------------------------------------

LaBGAScore_prep_s0_define_directories;


%% UNZIP IMAGES, SMOOTH, ZIP, SMOOTHED IMAGES, AND DELETE ALL UNZIPPED IMAGES
%----------------------------------------------------------------------------

if ~isempty(subjs2smooth)
    [C,ia,~] = intersect(derivsubjs,subjs2smooth);
    
    if ~isequal(C',subjs2smooth)
        error('subject defined in subjs2smooth not present in derivdir, please check');
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

cd(derivdir);