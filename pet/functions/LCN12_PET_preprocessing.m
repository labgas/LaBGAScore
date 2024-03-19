function go = LCN12_PET_preprocessing(dir_pet,infostring_tracer,infostring_PET,infostring_frames,dir_anat,name_logfile)
% assumptions:
%   Data-organisation:
%	    Data are transformed to nifti format and organized according to BIDS.
%       You have to specify a main directory where al the folders of the
%       subjects can be found. 
%       In the folder of each subject, there might be a folder ses-xx 
%       which contains the folder anat with the T1 weighted structural MRI 
%       and a folder pet which contains the dynamic PET images. 
%       If the folder ses-xx is not existing, we assume that the subfolder
%       anat and pet are directly under the main folder of the subject.
%       
%   PET acquisition, units, format, frame definition
%       A template for the PET naming should be specified. The PET data are 
%       named as 'subjectname*"infostring_PET"*.nii'.
%
%       PET data are acquired dynamically either from the start of the 
%       injection or as a number of frames between two time points.
%
%       The unit of the PET data is Bq/ml.
%
%       Images are decay corrected to the start of the injection/begin of scanning.
%
%       PET data are in 4D nifti format.
%
%       In the pet folder, there must be a .m file containing the frame
%       defintion by specifying the variable frames_timing which is a N x 2 
%       or N x 3 array (N = number of frames) for which the first column is 
%       the start time of the frame in seconds post injection, the second 
%       column is the end time of each frame in seconds post injection and 
%       the third column is optional with the weight for the frame (positive 
%       values). If weights are not specified, all frames are weighted in 
%       the same way. The frame definition file (in the same directory as 
%       the dynamic PET data) is named as 
%       subjectname_*"infostring_tracer"_"infostring_PET"_"infostring_frames".m
%
%       It has been checked that the CT or MRI (used for attenuation correction)
%       and the PET data are aligned.
%
%       After the previous check, the origin of the PET and MRI data have 
%       been set to AC and if necessary the pitch (and yaw and roll) have 
%       been adapted so that the initial orientation is more similar to the 
%       templates in MNI.
%       
% INPUT
%   dir_pet           = directory of the subject where PET data can be found
%   infostring_tracer = name used in the BIDS field trc (eg trc-NAV or trc-DPA714)
%   infostring_PET    = name used in the BIDS field rec (eg rec-dyn_PET)
%   infostring_frames = name used to indicate the frames (eg 'frames')
%   dir_anat          = directory where to find the T1 weighted structural MRI
%   OPTIONAL: name_logfile = full name of a logfile. If not specified, it
%             will be written to dir_pet and named 'LCN12_PET_preprocessing_log.txt' 
% Output:
%       in dir_pet the preprocessed data will be written, i.e.
%           warped dynamic data and summed images all in MNI space (having
%           a prefix 'w').
%       in the subdirectory tmp all intermediate data will be written
%       in the subdirectory batch all batch files will be written
%       in the subdirectory dir_anat the warped segmnentations of the MRI
%           coregistered to the PET are written
%
% author: Patrick Dupont, Lukas Van Oudenhove
% date: July 2023
% history: Nov 2023: we now use the LCN12_realign_dy_PET function for the
%                    realignment which works better for tracer which change
%                    their distribution too much over time. 
%          Mar 2024: made small edits to allow gzipped images to be used (LVO)
%                    added a line of bash code unannexing .nii files to prevent 
%                    writing permission issues on the LaBGAS server (LVO)
%
% THIS IS RESEARCH SOFTWARE
%__________________________________________________________________________
% @(#)LCN12_PET_preprocessing      v0.3           last modified: 2024/03/14
 
%------------- GENERAL SETTINGS -------------------------------------------
% limits for head movement
threshold_translation = 6; % maximum overall movement in mm
threshold_rotation    = 6; % maximum overall rotation in degrees
%++++ END OF SETTINGS - DO NOT CHANGE BELOW THIS LINE +++++++++++++++++++++

% determine the path of the prior data of SPM
tmp = which('spm.m');
[spm_pth,~,~]  = fileparts(tmp);
spm_pth_priors = fullfile(spm_pth,'tpm');
% find the name of the subject
ok = 0;
tmp1 = dir_pet;
while ok == 0
      [tmp1,name,~] = fileparts(tmp1);
      if strcmp(name(1:3),'sub') == 1
         subjectname = name;
         ok = 1;
      end
end
if exist('name_logfile','var') ~= 1 % logfile is not specified
   name_logfile = fullfile(dir_pet,'LCN12_PET_preprocessing_log.txt');      
end
fid  = fopen(name_logfile,'a+');
fprintf(fid,'subject = %s\n',subjectname);
fprintf(fid,'LCN12_PET_preprocessing.m \n');
fprintf(fid,'%c','-'*ones(1,30));
fprintf(fid,'\n');
fprintf(fid,'Processing started: %s \n',datetime('now'));
fprintf(fid,'Settings for maximum movement\n');
fprintf(fid,'\t threshold (overall translation in mm)   = %4.2f \n',threshold_translation);
fprintf(fid,'\t threshold (overall rotation in degrees) = %4.2f \n',threshold_rotation);
    
% define the directories we need
dir_tmp   = fullfile(dir_pet,'tmp');
dir_batch = fullfile(dir_pet,'batch');
dir_anat2 = fullfile(dir_pet,'anat');
    
% create the subdirectories batch and tmp in the folder dir_pet
cd(dir_pet)
mkdir(dir_tmp)
mkdir(dir_batch)
mkdir(dir_anat2) % contains the T1 data used for the PET processing
    
% check if all data are found for further processing
go = 1;
    
% read frame definition file
filelist = dir([subjectname '_*' infostring_tracer '_' infostring_PET '_' infostring_frames '.m']);
if size(filelist,1) ~= 1
   fprintf('ERROR subject %s: no or more than one frame defintion file found in %s\n',subjectname,dir_pet);
   fprintf(fid,'ERROR subject %s: no or more than one frame defintion file found in %s\n',subjectname,dir_pet);
   go = 0;
end
if go == 1
   frame_definition_file = filelist(1).name;
   fprintf(fid,'Timing file:  %s\n',frame_definition_file);
   copyfile(frame_definition_file,'tmp_frames.m');
   frames_timing = [];
   tmp_frames; % the variable frames_timing is now known
   delete('tmp_frames.m');
   nr_frames = size(frames_timing,1);
   % calculate midscan times (in minutes)
   acqtimes       = frames_timing(:,1:2)/60; % in minutes

   frame_5min = find(acqtimes(:,2) <= 5,1,'last');
   if frame_5min > 1 % assumed that the data are measured dynamically from the start of the injection
      go_dynamic_from_injection = 1;
      fprintf(fid,'assuming dynamic data measured from the start of the injection \n');
      fprintf(fid,'%i frames found starting at %i s after the injection until %i s post injection\n',nr_frames,frames_timing(1,1),frames_timing(end,2));
      fprintf('assuming dynamic data measured from the start of the injection \n');
      fprintf('%i frames found starting at %i s after the injection until %i s post injection\n',nr_frames,frames_timing(1,1),frames_timing(end,2));
   else % assumed that the data are measured dynamically between two later time points
      go_dynamic_from_injection = 0;        
      fprintf(fid,'assuming dynamic data measured between two time points\n');
      fprintf(fid,'%i frames found starting at %i s after the injection until %i s post injection\n',nr_frames,frames_timing(1,1),frames_timing(end,2));
      fprintf('assuming dynamic data measured between two time points \n');
      fprintf('%i frames found starting at %i s after the injection until %i s post injection\n',nr_frames,frames_timing(1,1),frames_timing(end,2));
   end

   cd(dir_pet)
   ! git annex unannex *.nii
   filelist = dir([subjectname '*' infostring_PET '*.nii']);
   if size(filelist,1) ~= 1
      fprintf('ERROR subject %s: expecting a single PET image (3D or 4D) in the directory PET \n',subjectname);
      fprintf(fid,'ERROR subject %s: expecting a single PET image (3D or 4D) in the directory PET \n',subjectname);
      go = 0;
   else
      Vref = spm_vol(filelist(1).name);
      if size(Vref,1) ~= nr_frames
         fprintf('ERROR subject %s: expecting %i frames but only %i found \n',subjectname,nr_frames,size(Vref,1));
         fprintf(fid,'ERROR subject %s: expecting %i frames but only %i found \n',subjectname,nr_frames,size(Vref,1));
         go = 0;
      else
         PET_filename = fullfile(dir_pet,filelist(1).name);
         fprintf(fid,'PET filename = %s \n',PET_filename); 
         fprintf(fid,'number of frames = %i \n',nr_frames);
      end
   end
end
    
% check if the structural T1w MRI is found
cd(dir_anat)
filelist = dir([subjectname '*_T1w.nii*']); % Lukas' edit to accommodate .nii.gz files
if size(filelist,1) ~= 1
   fprintf('ERROR subject %s: no or more than one structural scan found \n',subjectname);
   fprintf(fid,'ERROR subject %s: no or more than one structural scan found \n',subjectname);
   go = 0;
else
   filename  = filelist(1).name;
   T1w_image = fullfile(dir_anat2,filename);
   copyfile(fullfile(dir_anat,filename),T1w_image);
   if contains(T1w_image,'.gz')
       gunzip(T1w_image);
       delete(T1w_image);
       T1w_image = T1w_image(1,1:end-3);
   end
   fprintf(fid,'T1w filename = %s \n',T1w_image); 
end

if go ~= 1
   fprintf('not all data are found for subject %s. Skipping this subject \n',subjectname); 
   fprintf(fid,'not all data are found for subject %s. Skipping this subject \n',subjectname); 
else % all data are found
   %  initialize SPM
   fprintf('...initializing SPM \n');
   spm_jobman('initcfg')
   figure(1)
   % realign the PET data
   cd(dir_pet)
   fprintf(fid,'\n');
   fprintf(fid,'Realignment started: %s \n',datetime('now'));
   LCN12_realign_dy_PET(PET_filename,fullfile(dir_pet,frame_definition_file),threshold_translation,threshold_rotation);
   fprintf(fid,'Realignment ended at: %s\n',datetime('now'));

   % excluded frames
   [~,name,~] = fileparts(PET_filename);
   cd(dir_pet)
   tmp = load(['excluded_frames_' name '.mat'],'excluded_frames');
   excluded_frames = tmp(1).excluded_frames;
   fprintf('excluded frames PET based on excess head movement = ');
   fprintf(fid,'excluded frames PET based on excess head movement = ');
   for i = 1:length(excluded_frames)
       fprintf(fid,'%i ',excluded_frames(i)); 
       fprintf('%i ',excluded_frames(i)); 
   end
   fprintf(fid,'\n');
   fprintf('\n');
   fprintf('If excluded frames are present, evaluate if data are still reliable! \n');
   fprintf(fid,'If excluded frames are present, evaluate if data are still reliable! \n');
   
   % calculate a mean image
   if go_dynamic_from_injection ~= 1 % we will calculate the mean image excluding excluded frames
      filelist = cell(nr_frames,1);
      for i = 1:nr_frames
          filelist{i,1} = [PET_filename ',' num2str(i)];
      end        
      index_ok = setdiff(1:nr_frames,excluded_frames);
      filelist_ok = filelist(index_ok);
      matlabbatch = {};
      matlabbatch{1}.spm.spatial.realign.write.data = filelist_ok;
      matlabbatch{1}.spm.spatial.realign.write.roptions.which = [0 1];
      matlabbatch{1}.spm.spatial.realign.write.roptions.interp = 4;
      matlabbatch{1}.spm.spatial.realign.write.roptions.wrap = [0 0 0];
      matlabbatch{1}.spm.spatial.realign.write.roptions.mask = 1;
      matlabbatch{1}.spm.spatial.realign.write.roptions.prefix = 'r';
      % run batchfile
      spm_jobman('run',matlabbatch)
      filename_mean = fullfile(dir_pet,['mean' name '.nii']);
   end
   
   % coregister the PET and MRI
   %+++++++++++++++++++++++++++
   ref_image = T1w_image;
   if go_dynamic_from_injection ~= 1
      source_image = filename_mean;
   else
      filename_sum5min = fullfile(dir_tmp,'PETsum5min.nii');
      source_image = filename_sum5min;
   end
   filelist = cell(nr_frames,1);
   for i = 1:nr_frames
       filelist{i,1} = [PET_filename ',' num2str(i)];
   end        
   % genarate the batch file 
   fprintf('...creating the coregister batchfile for subject %s\n',subjectname);
   matlabbatch = {};
   matlabbatch{1}.spm.spatial.coreg.estimate.ref = cellstr(ref_image);
   matlabbatch{1}.spm.spatial.coreg.estimate.source = cellstr(source_image);
   matlabbatch{1}.spm.spatial.coreg.estimate.other = filelist;
   matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
   matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
   matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
   matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
   cd(dir_batch)
   save batch_coregister_MRI_PET matlabbatch
   fprintf(fid,'\n');
   fprintf(fid,'Coregistration started: %s \n',datetime('now'));
   % run batchfile
   spm_jobman('run',matlabbatch)
   fprintf(fid,'Coregistration ended at: %s\n',datetime('now'));

   cd(dir_anat2)
   % segment the MRI (which will also generate the warping to MNI)
   %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   matlabbatch = struct([]);
   matlabbatch{1}.spm.spatial.preproc.channel.vols = cellstr([T1w_image ',1']);
   matlabbatch{1}.spm.spatial.preproc.channel.biasreg = 0.001;
   matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 60;
   matlabbatch{1}.spm.spatial.preproc.channel.write = [0 1];
   matlabbatch{1}.spm.spatial.preproc.tissue(1).tpm = cellstr(fullfile(spm_pth_priors,'TPM.nii,1'));
   matlabbatch{1}.spm.spatial.preproc.tissue(1).ngaus = 1;
   matlabbatch{1}.spm.spatial.preproc.tissue(1).native = [1 0];
   matlabbatch{1}.spm.spatial.preproc.tissue(1).warped = [0 0];
   matlabbatch{1}.spm.spatial.preproc.tissue(2).tpm = cellstr(fullfile(spm_pth_priors,'TPM.nii,2'));
   matlabbatch{1}.spm.spatial.preproc.tissue(2).ngaus = 1;
   matlabbatch{1}.spm.spatial.preproc.tissue(2).native = [1 0];
   matlabbatch{1}.spm.spatial.preproc.tissue(2).warped = [0 0];
   matlabbatch{1}.spm.spatial.preproc.tissue(3).tpm = cellstr(fullfile(spm_pth_priors,'TPM.nii,3'));
   matlabbatch{1}.spm.spatial.preproc.tissue(3).ngaus = 2;
   matlabbatch{1}.spm.spatial.preproc.tissue(3).native = [1 0];
   matlabbatch{1}.spm.spatial.preproc.tissue(3).warped = [0 0];
   matlabbatch{1}.spm.spatial.preproc.tissue(4).tpm = cellstr(fullfile(spm_pth_priors,'TPM.nii,4'));
   matlabbatch{1}.spm.spatial.preproc.tissue(4).ngaus = 3;
   matlabbatch{1}.spm.spatial.preproc.tissue(4).native = [1 0];
   matlabbatch{1}.spm.spatial.preproc.tissue(4).warped = [0 0];
   matlabbatch{1}.spm.spatial.preproc.tissue(5).tpm = cellstr(fullfile(spm_pth_priors,'TPM.nii,5'));
   matlabbatch{1}.spm.spatial.preproc.tissue(5).ngaus = 4;
   matlabbatch{1}.spm.spatial.preproc.tissue(5).native = [1 0];
   matlabbatch{1}.spm.spatial.preproc.tissue(5).warped = [0 0];
   matlabbatch{1}.spm.spatial.preproc.tissue(6).tpm = cellstr(fullfile(spm_pth_priors,'TPM.nii,6'));
   matlabbatch{1}.spm.spatial.preproc.tissue(6).ngaus = 2;
   matlabbatch{1}.spm.spatial.preproc.tissue(6).native = [0 0];
   matlabbatch{1}.spm.spatial.preproc.tissue(6).warped = [0 0];
   matlabbatch{1}.spm.spatial.preproc.warp.mrf = 1;
   matlabbatch{1}.spm.spatial.preproc.warp.cleanup = 1;
   matlabbatch{1}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
   matlabbatch{1}.spm.spatial.preproc.warp.affreg = 'mni';
   matlabbatch{1}.spm.spatial.preproc.warp.fwhm = 0;
   matlabbatch{1}.spm.spatial.preproc.warp.samp = 3;
   matlabbatch{1}.spm.spatial.preproc.warp.write = [1 1];
   cd(dir_batch)
   save batch_segment matlabbatch
   cd(dir_anat2)
   fprintf(fid,'\n');
   fprintf(fid,'Segmentation started: %s \n',datetime('now'));
   % run batchfile
   spm_jobman('run',matlabbatch)
   fprintf(fid,'Segmentation ended at: %s\n',datetime('now'));
        
   % apply the warping parameters to the PET data
   %+++++++++++++++++++++++++++++++++++++++++++++
   [pth,filename,~] = fileparts(T1w_image);
   deformation_field = fullfile(pth,['y_' filename '.nii']);
   fprintf('...creating the warping batchfile for subject %s\n',subjectname);
   matlabbatch = {};
   matlabbatch{1}.spm.spatial.normalise.write.subj.def = cellstr(deformation_field);
   matlabbatch{1}.spm.spatial.normalise.write.subj.resample = filelist;
   matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70
                                                              78   76  85];
   matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = [2 2 2];
   matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 4;
   cd(dir_batch)
   save batch_warping_PET matlabbatch
   cd(dir_tmp)
   fprintf(fid,'\n');
   fprintf(fid,'Warping PET started: %s \n',datetime('now'));
   % run batchfile
   spm_jobman('run',matlabbatch)
   fprintf(fid,'Warping PET ended at: %s\n',datetime('now'));
    
   % apply the warping parameters to the GM and WM images
   %+++++++++++++++++++++++++++++++++++++++++++++++++++++   
   fprintf('...creating the warping batchfile for subject %s\n',subjectname);
   filelistMRI = {};
   filelistMRI{1,1} = fullfile(pth,['c1' filename '.nii,1']); % GM
   filelistMRI{2,1} = fullfile(pth,['c2' filename '.nii,1']); % WM
   filelistMRI{3,1} = fullfile(pth,['c3' filename '.nii,1']); % CSF
   filelistMRI{4,1} = fullfile(pth,['m' filename '.nii,1']); % bias corrected MRI
   matlabbatch = {};
   matlabbatch{1}.spm.spatial.normalise.write.subj.def = cellstr(deformation_field);
   matlabbatch{1}.spm.spatial.normalise.write.subj.resample = filelistMRI;
   matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70
                                                              78   76  85];
   matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = [1 1 1];
   matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 4;
   cd(dir_batch)
   save batch_warping_GM_WM_CSF matlabbatch
   cd(dir_tmp)
   fprintf(fid,'\n');
   fprintf(fid,'Warping segmentations started: %s \n',datetime('now'));
   % run batchfile
   spm_jobman('run',matlabbatch)
   fprintf(fid,'Warping segmentations ended at: %s\n',datetime('now'));
end
fclose(fid);
end