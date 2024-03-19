function LCN12_realign_dy_PET(filename_dyPET,frame_definition_file,threshold_translation,threshold_rotation)
% This function will realign dynamic PET images. 
% If measured from the start of the injection, the first 5 minutes will be 
% summed and it is assumed that there is not movement during this 5 minutes 
% (which has to be checked) and that the attenuation correction file is 
% alligned with this image (this can be checked by using the CT or MRI upon 
% which the attenuation correcton file is based).
% If not measured from the start of the injection, we assume that the first
% dynamic frame is alligned with the attenuation correcton file). 
%
% INPUT
%   filename_dyPET = full file name of the dynamic PET data (assumed that
%                    they are in 4D nifti)
%   frame_definition_file = full file name of the .m file in which the 
%                    variable frames_timing is defind which is a n x 3 matrix 
%                    with n the number of frames and the first column is the
%                    start time of the frame; the second column is the end 
%                    time of the frame and the third column (optional) is 
%                    the weight for the frame. All times are in seconds 
%                    post-injection
%   threshold_translation = maximum translation in x, y or z in mm
%   threshold_rotation    = maximum rotation in roll, yaw and pitch in degrees
% 
% OUTPUT
%   the transformation in the headers of the dynamic PET file will be
%       adapted (so make a copy first if you want to preserve the 
%       "untouched" data)
%   a sum5min.nii image will be generated if measured from the start of the
%       injection.
%   a mat file containinging the bad_frames will be saved. IMPORTANT: this
%   is corresponding to the effective frame numbers (excluded_frames.mat).
%   
% author: Patrick Dupont
% date: June 2021
% history: June 2021: a file excluded_frames.mat is now saved which contains
%                     the effective frame numbers to be excluded later on.
%          Sept 2021: we changed the coregistration algorithm from
%                     normalized mutual information to normalized cross 
%                     correlation. The latter works much better for data 
%                     such as PET MK6240.
%          Nov 2023: BIDS compatible
%
% THIS IS RESEARCH SOFTWARE
%__________________________________________________________________________
% @(#)LCN12_realign_dy_PET.m        v0.3          last modified: 2023/11/09

% find the directories where to find the specified files
[dir_pet,name,~] = fileparts(filename_dyPET);

dir_tmp   = fullfile(dir_pet,'tmp');
dir_batch = fullfile(dir_pet,'batch');

% read frame definition file
if ~isempty(dir_pet)
   cd(dir_pet)
end
% matlab cannot handdle - in the name so we copy it to a tmp file
copyfile(frame_definition_file,'frames.m');
frames_timing = [];
frames; % the variable frames_timing is now known
delete('frames.m');

nr_frames = size(frames_timing,1);
% calculate midscan times (in minutes)
acqtimes       = frames_timing(:,1:2)/60; % in minutes
frame_5min = find(acqtimes(:,2) <= 5,1,'last');

matlabbatch = {};
spm_jobman('initcfg')
if ~isempty(dir_pet)
   cd(dir_pet)
end

Vref = spm_vol(filename_dyPET);
nr_frames_dy = size(Vref,1);
Pdy = cell(nr_frames_dy,1);
if nr_frames_dy ~= size(frames_timing,1)
   fprintf('ERROR: total number of frames is not consistent with the frame defintion file\n');
   return;
else
   for frame = 1:nr_frames_dy
       Pdy{frame} = [filename_dyPET ',' num2str(frame)];
   end
end    

if frame_5min > 1 % assumed that the data are measured dynamically from the start of the injection
   go_dynamic_from_injection = 1;
   
   % make sum image of the first 5 min 
   %++++++++++++++++++++++++++++++++++
   img5min = zeros(Vref(1).dim(1),Vref(1).dim(2),Vref(1).dim(3));
   for i = 1:frame_5min
       [img,~] = LCN12_read_image(char(Pdy(i)),Vref(1));
       img5min = img5min + img;
   end
   filename_sum5min = fullfile(dir_tmp,'PETsum5min.nii');
   LCN12_write_image(img5min,filename_sum5min,'sum first 5 min',Vref(1).dt(1),Vref(1));        
else % assumed that the data are measured dynamically between two later time points
   go_dynamic_from_injection = 0;        
end

if go_dynamic_from_injection == 1
   fid = fopen('rp_sum5min.txt','w+');
   Pdy2 = Pdy;
   Pdy2{frame_5min} = [filename_sum5min ',1'];
   % perform coregistration with previous frame
   tmp0 = spm_imatrix(Vref(frame_5min).mat);
   for i = frame_5min+1:nr_frames_dy
       clear ref_image source_image Vorig Vnew tmp_orig tmp_new
       ref_image = char(Pdy2{i-1});
       source_image = char(Pdy2{i});
       % generate the batch file 
       matlabbatch = {};
       matlabbatch{1}.spm.spatial.coreg.estimate.ref = cellstr(ref_image);
       matlabbatch{1}.spm.spatial.coreg.estimate.source = cellstr(source_image);
%        matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
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
   save batch_realign_last_pair matlabbatch
   fclose(fid);
elseif go_dynamic_from_injection == 0
   filelist = cell(nr_frames_dy,1);
   for i = 1:nr_frames_dy
       filelist{i,1} = char(Pdy{i});
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
if go_dynamic_from_injection == 1
   excluded_frames = bad_frames + frame_5min;
   eval(['save excluded_frames_' name ' excluded_frames ']);   
elseif go_dynamic_from_injection == 0
   excluded_frames = bad_frames;
   eval(['save excluded_frames_' name ' excluded_frames ']);
end
% read all the frames to generate an image for quality checking
img4D = zeros(Vref(1).dim(1),Vref(1).dim(2),Vref(1).dim(3),nr_frames);
for i = 1:nr_frames_dy
    img = LCN12_read_image(char(Pdy(i)),Vref(1));
    img4D(:,:,:,i) = img;
end

% generate images for quality checking
LCN12_generate_qc_from4D(img4D,dir_tmp);
end