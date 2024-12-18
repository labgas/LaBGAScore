function LCN12_realign2_dy_PET(PET_masked,PET_filename,filename_sum5min,frame_5min,threshold_translation,threshold_rotation)
% This function will realign dynamic PET images by realignment of masked 
% dynamic PET and it will be applied to the original dynamic PET images.  
% If measured from the start of the injection, a sum image of the first 5 
% minutes will be used and it is assumed that there is not movement during 
% this 5 minutes (which has to be checked) and that the attenuation 
% correction file is alligned with the first part of the dynamic PET data 
% (which has to be checked before preprocessing of the PET data starts).
% If not measured from the start of the injection, we assume that the first
% dynamic frame is alligned with the attenuation correcton file and 
% filename_sum5min should be '' and frame_5min should be 0. 
%
% INPUT 
%   PET_masked = full file name of the masked dynamic PET data (assumed that
%                    they are in 4D nifti)
%   PET_filename = full file name of the dynamic PET data (assumed that
%                    they are in 4D nifti)
%   filename_sum5min = full file name of the sum5min image (for dynamic
%                      data starting later after injection, give this the 
%                      value '').
%   frame_5min = last frame included in the sum5min image (for dynamic
%                      data starting later after injection, give this the 
%                      value 0). 
%   threshold_translation = maximum translation in x, y or z in mm
%   threshold_rotation    = maximum rotation in roll, yaw and pitch in degrees
% 
% OUTPUT
%   the transformation in the headers of the dynamic PET files PET_masked and PET_filename
%       will be adapted (so make a copy first if you want to preserve the 
%       "untouched" data)
%   a mat file containinging the bad_frames will be saved. IMPORTANT: this
%   is corresponding to the effective frame numbers (excluded_frames.mat).
%   
% author: Patrick Dupont
% date: December 2024
% history:
%
% THIS IS RESEARCH SOFTWARE
%__________________________________________________________________________
% @(#)LCN12_realign2_dy_PET.m       v0.1          last modified: 2024/12/12

% find the directories where to find the specified files
[dir_pet,name,~] = fileparts(PET_filename);

dir_tmp   = fullfile(dir_pet,'tmp');
dir_batch = fullfile(dir_pet,'batch');

Vref = spm_vol(PET_masked);
nr_frames_dy = size(Vref,1);
Pdy_A  = cell(nr_frames_dy,1);
Pdy_B = cell(nr_frames_dy,1);
for frame = 1:nr_frames_dy
    Pdy_A{frame} = [PET_masked ',' num2str(frame)];
    Pdy_B{frame} = [PET_filename ',' num2str(frame)];
end
matlabbatch = {};
spm_jobman('initcfg')
if ~isempty(dir_pet)
   cd(dir_pet)
end
% read a dynamic dataset before realignment and generate image for QC
img4D = zeros(Vref(1).dim(1),Vref(1).dim(2),Vref(1).dim(3),nr_frames_dy);
for i = 1:nr_frames_dy
    img = LCN12_read_image(char(Pdy_A(i)),Vref(1));
    img4D(:,:,:,i) = img;
end
LCN12_generate_qc_from4D(img4D,dir_tmp,'qc_before_realignment_');

% realign data
if frame_5min > 1 % assumed that the data are measured dynamically from the start of the injection
   fid = fopen('rp_sum5min.txt','w+');
   Pdy2 = Pdy_A;
   Pdy2{frame_5min} = [filename_sum5min ',1'];
   % perform coregistration with previous frame
   tmp0 = spm_imatrix(Vref(frame_5min).mat);
   for i = frame_5min+1:nr_frames_dy
       clear ref_image source_image Vorig Vnew tmp_orig tmp_new
       ref_image = char(Pdy2{i-1});
       source_image = char(Pdy2{i});
       other_image  = char(Pdy_B{i});
       % generate the batch file 
       matlabbatch = {};
       matlabbatch{1}.spm.spatial.coreg.estimate.ref    = cellstr(ref_image);
       matlabbatch{1}.spm.spatial.coreg.estimate.source = cellstr(source_image);
       matlabbatch{1}.spm.spatial.coreg.estimate.other  = cellstr(other_image);
       matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
       matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'ncc';
       matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
       matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
       matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
       % run batchfile
       spm_jobman('run',matlabbatch)
       Vnew = spm_vol(source_image);
       tmp_new  = spm_imatrix(Vnew.mat);
       fprintf(fid,'%8.6f \t %8.6f \t %8.6f \t %8.6f \t %8.6f \t %8.6f \n',tmp_new(1:6) - tmp0(1:6));
   end
   cd(dir_batch);
   save batch_realign_pairwise matlabbatch
   fclose(fid);
elseif frame_5min == 0
   filelist = cell(nr_frames_dy,1);
   for i = 1:nr_frames_dy
       filelist{i,1} = char(Pdy_A{i});
   end
   % generate the batch file for the realignment
   matlabbatch = {};
   matlabbatch{1}.spm.spatial.realign.estimate.data = {filelist}';
   matlabbatch{1}.spm.spatial.realign.estimate.eoptions.quality = 0.9;
   matlabbatch{1}.spm.spatial.realign.estimate.eoptions.sep = 4;
   matlabbatch{1}.spm.spatial.realign.estimate.eoptions.fwhm = 5;
   matlabbatch{1}.spm.spatial.realign.estimate.eoptions.rtm = 1;
   matlabbatch{1}.spm.spatial.realign.estimate.eoptions.interp = 2;
   matlabbatch{1}.spm.spatial.realign.estimate.eoptions.wrap = [0 0 0];
   matlabbatch{1}.spm.spatial.realign.estimate.eoptions.weight = '';
   cd(dir_batch);
   save batch_realign matlabbatch
   % run batchfile
   spm_jobman('run',matlabbatch)
end

% display realignment parameters
cd(dir_pet)
filelist = dir('rp_*.txt');
if size(filelist,1) ~= 1
   fprintf('ERROR subject %s: more than 1 realignment file (rp_*.txt) found\n',name);
end
rp_file = filelist(1).name;
bad_frames = LCN12_analyze_headmovement_PET(rp_file,threshold_translation,threshold_rotation);

% transform bad_frames to excluded_frames containing the effective frame
% numbers
if frame_5min > 1
   excluded_frames = bad_frames + frame_5min;
   eval(['save excluded_frames_' name ' excluded_frames ']);   
elseif frame_5min == 0
   excluded_frames = bad_frames;
   eval(['save excluded_frames_' name ' excluded_frames ']);
end
% read all the frames to generate an image for quality checking
img4D = zeros(Vref(1).dim(1),Vref(1).dim(2),Vref(1).dim(3),nr_frames_dy);
for i = 1:nr_frames_dy
    img = LCN12_read_image(char(Pdy_A(i)),Vref(1));
    img4D(:,:,:,i) = img;
end
LCN12_generate_qc_from4D(img4D,dir_tmp,'qc_after_realignment_');
end