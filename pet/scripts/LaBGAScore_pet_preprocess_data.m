% discoverie_pet_s1_preprocess_data.m
%
% assumptions:
%   Data-organisation:
%	    Data are transformed to nifti format and organized according to BIDS.
%       You have to specify a main directory where al the folders of the
%       subjects can be found. 
%       In the folder of each subject, there might be a folder ses-xx 
%       which contains the folder anat with the T1 weighted structural MRI 
%       and a folder pet_"infostring_tracer" which contains the PET data 
%       measured with the specified tracer. If the folder ses-xx is not 
%       existing, we assume that the subfolder anat and 
%       pet_"infostring_tracer" are directly under the main folder of the subject.
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
%               
% author: Patrick Dupont, Lukas Van Oudenhove
% date: September 2023
% history: October 2023 cleaning code
%                       taking into account the naming of the pet folder
%          March 2024   adapted to work with standard LaBGAS file
%                       organization (LVO)
%                       adapted to move output to derivatives
%
% THIS IS RESEARCH SOFTWARE
%__________________________________________________________________________
% @(#)LCN12_PET_preprocess_data.m       v0.13     last modified: 2024/03/15

clear
close all

LaBGAScore_prep_s0_define_directories; % change to study-specific script

maindir           = BIDSdir; % directory where the folders of each subject can  be found
sessiondir        = ''; % if empty, we assume that there is no folder session and the folders anat and pet are directly under the subject folder
infostring_tracer = 'trc-DPA714';
infostring_PET    = 'rec-acdyn_pet'; % the PET data are "subjectname"_"sessiondir"_"infostring_tracer"_"infostring_PET".nii (example =  sub-test_trc-DPA714_rec-acdyn_pet.nii)
infostring_frames = 'frames'; % in the folder pet, we assume a .m file name "subjectname"_"sessiondir"_"infostring_tracer"_"infostring_PET"_"infostring_frames".m (example sub-test_trc-DPA714_rec-acdyn_pet_frames.m)


% Lukas' code to automate selection of subjects who have a pet dir in BIDS
% USE THIS FOR PREPROCESSING ALL PET SUBJECTS AT ONCE
dir_BIDS = dir(fullfile(BIDSdir,'sub-*'));

petsubcounter = 1;

for sub = 1:size(dir_BIDS,1)
    BIDSsubdir = fullfile(BIDSdir,dir_BIDS(sub).name);
    dir_BIDSsubdir = dir(BIDSsubdir);
    if contains([dir_BIDSsubdir(:).name],'pet')
        SUBJECTS{petsubcounter,1} = dir_BIDS(sub).name;
        petsubcounter = petsubcounter + 1;
    else
        continue
    end
end

% Patrick's code to enter subjects manually
% USE THIS FOR PREPROCESSING SELECTED SUBJECTS
% SUBJECTS = {
% subjectdir
% 'sub-KUL048'
% 'sub-KUL059'
% 'sub-KUL102'
% 'sub-KUL102B'
% 'sub-KUL117'
% };

% maindir           = 'C:\DATA\LAURE\NAV'; % directory where the folders of each subject can  be found
% sessiondir        = 'ses-01'; % if empty, we assume that there is no folder session and the folders anat and pet are directly under the subject folder
% infostring_tracer = 'trc-NAV4694';
% infostring_PET    = 'rec-dyn_pet'; % the PET data are "subjectname"_"sessiondir"_"infostring_tracer"_"infostring_PET".nii (example =  sub-900_ses-01_trc-NAV_rec-dyn_pet.nii)
% infostring_frames = 'frames'; % in the folder pet, we assume a .m file name "subjectname"_"sessiondir"_"infostring_tracer"_"infostring_PET"_"infostring_frames".m (example NAV_frames.nii)
% 
% SUBJECTS = {
% % subjectdir
% 'sub-894'
% 'sub-900'
% % 'sub-901'
% };

% maindir           = 'C:\DATA\LAURE\MK6240'; % directory where the folders of each subject can  be found
% sessiondir        = 'ses-01'; % if empty, we assume that there is no folder session and the folders anat and pet are directly under the subject folder
% infostring_tracer = 'trc-MK6240';
% infostring_PET    = 'rec-dyn_pet'; % the PET data are "subjectname"_"sessiondir"_"infostring_tracer"_"infostring_PET".nii (example =  sub-900_ses-01_trc-NAV_rec-dyn_pet.nii)
% infostring_frames = 'frames'; % in the folder pet, we assume a .m file name "subjectname"_"sessiondir"_"infostring_tracer"_"infostring_PET"_"infostring_frames".m (example NAV_frames.nii)
% 
% SUBJECTS = {
% % subjectdir
% 'sub-794'
% 'sub-798'
% };

%++++ END OF SETTINGS - DO NOT CHANGE BELOW THIS LINE +++++++++++++++++++++

nr_subjects = size(SUBJECTS,1);

% determine the path of the prior data of SPM
tmp             = which('spm.m');
[spm_pth,~,~]   = fileparts(tmp);
spm_pth_priors  = fullfile(spm_pth,'tpm');

curdir = pwd;
for subj = 1:nr_subjects
    clear subjectdir subjectname dir_anat dir_pet name_logfile fid go0   
   
    subjectname = SUBJECTS{subj};
    subjectdir  = fullfile(maindir,subjectname);
    if ~isempty(sessiondir)
       dir_anat  = fullfile(fullfile(subjectdir,sessiondir),'anat'); 
       dir_pet   = fullfile(fullfile(subjectdir,sessiondir),['pet_' infostring_tracer]);
    else
       dir_anat = fullfile(subjectdir,'anat'); 
       dir_pet  = fullfile(subjectdir,['pet_' infostring_tracer]); 
    end
    
    fprintf('working on subject %s \n',subjectname);
    name_logfile = fullfile(dir_pet,'LCN12_PET_preprocess_data_log.txt');
    fid  = fopen(name_logfile,'a+');
    fprintf(fid,'subject = %s\n',subjectname);
    fprintf(fid,'LCN12_PET_preprocess_data.m \n');
    fprintf(fid,'%c','-'*ones(1,30));
    fprintf(fid,'\n');
    fprintf(fid,'Processing started at %s\n',datetime('now'));
    go0 = LCN12_PET_preprocessing(dir_pet,infostring_tracer,infostring_PET,infostring_frames,dir_anat,name_logfile);
    if go0 == 1
       fprintf(fid,'Processing ended at %s\n',datetime('now'));
    else
       fprintf(fid,'problem with preprocessing subject %s. Please check! \n',subjectname);
    end
    fprintf(fid,'++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n\n');
    fclose(fid);
    
    % move to derivatives - Lukas' code
    list_pet = dir(dir_pet);
    list_pet2copy = list_pet(~startsWith({list_pet(:).name},subjectname));
    derivrootdir = fileparts(derivdir);
    if ~isempty(sessiondir)
        subjectderivpetdir = fullfile(derivrootdir,['pet_' infostring_tracer],subjectname,sessiondir);
    else
        subjectderivpetdir = fullfile(derivrootdir,['pet_' infostring_tracer],subjectname);
    end
    if ~exist(subjectderivpetdir,'dir')
        mkdir(subjectderivpetdir);
    end
    cd(dir_pet);
    for i = 3:size(list_pet2copy,1)
        movefile(fullfile(list_pet2copy(i).folder,list_pet2copy(i).name),subjectderivpetdir);
    end
    matfile2copy = dir(fullfile(list_pet2copy(i).folder,'*.mat'));
    movefile(fullfile(matfile2copy.folder,matfile2copy.name),subjectderivpetdir);
end

