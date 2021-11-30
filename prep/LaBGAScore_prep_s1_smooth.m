%% WAD_prep_s1_deriv_unzip_nii_smooth
%
% This script will unzip fMRIprep output images, smooth them, zip the
% smoothed images, and delete all the unzipped images again
% 
% DEPENDENCIES
% SPM12 on your Matlab path
% 
% INPUTS
% preprocessed .nii.gz images outputted by fMRIprep
%
% OUTPUT
% smoothed .nii.gz images
%
%__________________________________________________________________________
%
% author: Lukas Van Oudenhove and Iris Coppieters
% date:   March, 2021
%
%__________________________________________________________________________
% @(#)% WAD_prep_s1_deriv_unzip_nii_smooth.m         v1.1       
% last modified: 2021/06/04
%
% changes versus version 1.0
% built in option to specify a list of subjects to smooth, in addition to
% looping over all subjects

%% DEFINE DIRECTORIES, SMOOTHING OPTIONS, AND SUBJECTS
%--------------------------------------------------------------------------
derivdir = 'C:\Users\lukas\Dropbox (Dartmouth College)\4.Iris_WAD_MRI\derivatives\fmriprep'; %directory with fMRIprep output
fwhm = 6; % kernel width in mm
prefix = 's6-'; % prefix for name of smoothed images

subjs2smooth = {'sub-P010','sub-P012','sub-P013','sub-P016','sub-P018','sub-P019','sub-P023','sub-P027','sub-P028','sub-P029','sub-P050','sub-P070','sub-P106'}; % specify if you only want to smooth a subset of all subjects in derivdir, otherwise leave cell array empty

subjs=dir(fullfile(derivdir,'sub-*'));
subjdirs=char(subjs(:).name);
idx=[subjs(:).isdir]';
subjdirs=deblank(subjdirs(idx,:));

% DO NOT CHANGE CODE BELOW THIS LINE
% ALWAYS MAKE A LOCAL COPY OF EXAMPLE SCRIPTS BEFORE MODIFYING


%% UNZIP IMAGES, SMOOTH, ZIP, SMOOTHED IMAGES, AND DELETE ALL UNZIPPED IMAGES
%----------------------------------------------------------------------------
if ~isempty(subjs2smooth)
    [C,ia,~] = intersect(subjdirs,subjs2smooth);
    if ~isequal(C',subjs2smooth)
        error('subject defined in subjs2smooth not present in derivdir, please check');
    else
        for i=ia'
            cd(fullfile(derivdir,subjdirs(i,:),'/ses-1/func'));
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
            eval(['save ' subjdirs(i,:) '_smooth.mat matlabbatch']); 
            spm_jobman('initcfg');
            spm_jobman('run',matlabbatch);
            % zip smoothed files
            gzip('s6*');
            % delete all unzipped files
            delete('*.nii');
        end % for loop over subjs2smooth
    end % if loop checking intersection of subjs2smooth and subjdirs
else
    for i=1:size(subjdirs,1)
        cd(fullfile(derivdir,subjdirs(i,:),'/ses-1/func'));
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
        eval(['save ' subjdirs(i,:) '_smooth.mat matlabbatch']); 
        spm_jobman('initcfg');
        spm_jobman('run',matlabbatch);
        % zip smoothed files
        gzip('s6*');
        % delete all unzipped files
        delete('*.nii');
    end % for loop over subjdirs
end % if loop checking smoothing option