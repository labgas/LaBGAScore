% LaBGAScore_pet_s2_model_TSPO_DPA714.m
%
% assumptions: 
%	    Data are organized according to BIDS.
%       You have to specify a main directory where al the folders of the
%       subjects can be found. 
%       In the folder of each subject, there might be a folder ses-xx 
%       which contains the folder anat with the T1 weighted structural MRI 
%       and a folder pet which contains the dynamic PET images. 
%       If the folder ses-xx is not existing, we assume that the subfolder
%       anat and pet are directly under the main folder of the subject.
%       
%       The data are preprocessed using LCN12_PET_preprocess_data.m
%
%       A template for the PET naming should be specified. PET data are 
%       acquired dynamically from the start of the injection. The unit of 
%       the PET data is expresses as Bq/ml.Images are decay corrected to 
%       the start of the injection/begin of scanning.
%       We assume that PET data are in 4D nifti format.
%       In the pet folder, there must be a .m file containing the frame
%       defintion by specifying the variable frames_timing which is a N x 2 
%       or N x 3 array (N = number of frames) for which the first column is 
%       the start time of the frame in seconds post injection, the second 
%       column is the end time of each frame and the third column is 
%       optional with the weight for the frame (positive values). If
%       weights are not specified, all frames are weighted in the same way.
%               
% author: Patrick Dupont, Lukas Van Oudenhove
% date: October 2023
% history: february 2024: if the region based analysis is selected, an
%                         excel file is written which contains for all 
%                         subjects and all regions, the volume of distribution
%                         and the error of the fit.
%          March 2024:    adapted to work with standard LaBGAS file
%                         organization (LVO) and adapted version of
%                         preprocessing script
%                         LaBGAScore_pet_preprocess_data
%                         built in check for enough frames after
%                         logan_start_time (LVO)
%          May 2024       further adaptations to fit LaBGAS file organization (LVO)
%                         built in more options for automatic subject and ROI definition (LVO)
%                         fixed bug in writing of summary excel table for VOI results (LVO)
%
% THIS IS RESEARCH SOFTWARE
%__________________________________________________________________________
% @(#)LCN12_PET_TSPO_DPA714.m        v0.5         last modified: 2024/05/31

clear
close all

LaBGAScore_prep_s0_define_directories;

maindir           = BIDSdir; % directory where the folders of each subject can  be found
sessiondir        = ''; % if empty, we assume that there is no folder session and the folders anat and pet are directly under the subject folder
infostring_tracer = 'trc-DPA714';
infostring_PET    = 'rec-acdyn_pet'; % the PET data are "subjectname"_"sessiondir"_"infostring_tracer"_"infostring_PET".nii (example =  sub-test_trc-DPA714_rec-acdyn_pet.nii)
infostring_frames = 'frames'; % in the folder pet, we assume a .m file name "subjectname"_"sessiondir"_"infostring_tracer"_"infostring_PET"_"infostring_frames".m (example sub-test_trc-DPA714_rec-acdyn_pet_frames.m)
% infostring_input and infostring_metab need only be defined if we use
% models which require an arterial input function
infostring_input = 'data_blood'; % % in the folder pet, we assume a .m file name "subjectname"_"sessiondir"_"infostring_tracer"_"infostring_input".m (example sub-test_trc-DPA714_rec-acdyn_data_blood.nii)
infostring_metab = 'data_metab'; % % in the folder pet, we assume a .m file name "subjectname"_"sessiondir"_"infostring_tracer"_"infostring_metab".m (example sub-test_trc-DPA714_rec-acdyn_data_metab.nii)

% Lukas' code to take into account that preprocessed data are in
% derivatives/pet-<infostring_tracer>
derivrootdir = fileparts(derivdir);
derivpetdir = fullfile(derivrootdir,['pet_' infostring_tracer]);
firstlevelpetdir = fullfile(rootdir,'firstlevel',['pet_' infostring_tracer]);
    if ~exist(firstlevelpetdir,'dir')
        mkdir(firstlevelpetdir);
    end
secondlevelpetdir = fullfile(rootdir,'secondlevel',['pet_' infostring_tracer]);
    if ~exist(secondlevelpetdir,'dir')
        mkdir(secondlevelpetdir);
    end
secondlevelpetmaskdir = fullfile(secondlevelpetdir,'masks');
    if ~exist(secondlevelpetmaskdir,'dir')
        mkdir(secondlevelpetmaskdir);
    end   
secondlevelpetresultsdir = fullfile(secondlevelpetdir,'results');
    if ~exist(secondlevelpetresultsdir,'dir')
        mkdir(secondlevelpetresultsdir);
    end  
secondlevelpetroidir = fullfile(secondlevelpetmaskdir,'rois');
    if ~exist(secondlevelpetroidir,'dir')
        mkdir(secondlevelpetroidir);
    end 

% Lukas' code to automate selection of subjects who have a pet dir in BIDS
% USE THIS FOR ANALYZING ALL PET SUBJECTS AT ONCE
% dir_BIDS = dir(fullfile(BIDSdir,'sub-*'));
% 
% petsubcounter = 1;
% 
% for sub = 1:size(dir_BIDS,1)
%     BIDSsubdir = fullfile(BIDSdir,dir_BIDS(sub).name);
%     dir_BIDSsubdir = dir(BIDSsubdir);
%     if contains([dir_BIDSsubdir(:).name],'pet')
%         SUBJECTS{petsubcounter,1} = dir_BIDS(sub).name;
%         petsubcounter = petsubcounter + 1;
%     else
%         continue
%     end
% end
% 
% dir_derivpet = dir(fullfile(derivpetdir,'sub-*'));
% 
% if ~isequal(SUBJECTS,{dir_derivpet(:).name}')
%     error('#subjects with PET data in BIDS and derivatives subdatasets does not match, please preprocess all subjects before proceeding');
% end;

% Lukas' code to automate selection of subjects in derivatives/pet
% USE THIS FOR ANALYZING ALL PET SUBJECTS AT ONCE WITHOUT CHECKING CONSISTENCY BETWEEN BIDS AND DERIVATIVES SUBDATASETS
dir_derivpet = dir(fullfile(derivpetdir,'sub-*'));
SUBJECTS = {dir_derivpet(:).name}';

% Patrick's code to enter subjects manually
% USE THIS FOR ANALYZING SELECTED SUBJECTS
% SUBJECTS = {
% subjectdir
% 'sub-KUL048'
% 'sub-KUL061'
% 'sub-KUL102'
% 'sub-KUL102B'
% 'sub-KUL117'
% };

outputfile_excel_region_based = fullfile(secondlevelpetresultsdir,['results_' infostring_tracer '_VOIS.xlsx']);

go_input = 1; % 1 if we require an input function, 0 else
region_based_analysis = 1; % if 1, we expect that  you have specified a number of regions and the analysis will be performed in these regions
voxel_based_analysis  = 1; % if 1, we perform a voxel based modelling
figures_on   = 1; % if 1 and if region_based_analysis = 1, we show the figures for each region and pause. You need to hit a key to proceed.
save_figures = 1; % if 1 and if region_based_analysis = 1, we save the figures to file and pausing of figures is overruled.
additional_smooth_parametric = 8; % kernel size in mm; isotropic Gaussian 3D smoothing; additional smoothing (if a voxel based analysis is needed)

% list of VOIS in MNI space; 

% Lukas' code to select all ROIs from secondlevel/pet/masks/rois

dir_secondlevelpetroidir = dir(fullfile(secondlevelpetroidir,'*.nii'));

VOIS = cell(size(dir_secondlevelpetroidir,1),2);

for roi = 1:size(VOIS,1)
    VOIS{roi,1} = fullfile(dir_secondlevelpetroidir(roi).folder,dir_secondlevelpetroidir(roi).name);
    VOIS{roi,2} = dir_secondlevelpetroidir(roi).name(1:end-4);
end

% Patrick's code to specify paths to ROIs manually

% VOIS = { 
%     % full name (extensions .img or .nii)            short name of the VOI
%     'C:\DATA\ATLAS\VOIS_preclinAD\resl_to_wmc1_SPM12_aal_ant_cingulate.img'        'ant cingulate'
%     'C:\DATA\ATLAS\VOIS_preclinAD\resl_to_wmc1_SPM12_aal_composite_cortical.img'   'composite cortical'
%     'C:\DATA\ATLAS\VOIS_preclinAD\resl_to_wmc1_SPM12_aal_frontal.img'              'frontal'
%     'C:\DATA\ATLAS\VOIS_preclinAD\resl_to_wmc1_SPM12_aal_lat_temporal.img'         'lat temporal'
%     'C:\DATA\ATLAS\VOIS_preclinAD\resl_to_wmc1_SPM12_aal_med_temporal.img'         'med temporal'
%     'C:\DATA\ATLAS\VOIS_preclinAD\resl_to_wmc1_SPM12_aal_nucl_caudatus.img'        'nucl caudatus'
%     'C:\DATA\ATLAS\VOIS_preclinAD\resl_to_wmc1_SPM12_aal_occipital.img'            'occipital'     
%     'C:\DATA\ATLAS\VOIS_preclinAD\resl_to_wmc1_SPM12_aal_parietal.img'             'parietal'
%     'C:\DATA\ATLAS\VOIS_preclinAD\resl_to_wmc1_SPM12_aal_post_cingulate.img'       'post cingulate'
%     'C:\DATA\ATLAS\VOIS_preclinAD\resl_to_wmc1_SPM12_aal_putamen.img'              'putamen'
%     'C:\DATA\ATLAS\VOIS_preclinAD\resl_to_wmc1_SPM12_aal_striatum.img'             'striatum'
%     'C:\DATA\ATLAS\VOIS_preclinAD\resl_to_wmc1_SPM12_aal_thalamus.img'             'thalamus'
%     'C:\DATA\ATLAS\VOIS_preclinAD\resl_to_wmc1_SPM12_aal_whole_cerebellum.img'     'whole cerebellum'
%     };

%------------- GENERAL SETTINGS -------------------------------------------
logan_start_time   = 31; % in min (ref Van Weehaeghe et al. J Nucl Med. 2020 Apr;61(4):604-607)
p0_hill        = [50 -1];      % see LCN_calc_intact_tracer_hill for details

%++++ END OF SETTINGS - DO NOT CHANGE BELOW THIS LINE +++++++++++++++++++++

global TIME_METAB
global FRACTION_INTACT_TRACER
global WEIGHTS_METAB

nr_subjects = size(SUBJECTS,1);
nr_VOIS = size(VOIS,1);

results_region_based_VD    = cell(nr_subjects+1,nr_VOIS+1);
results_region_based_error = cell(nr_subjects+1,nr_VOIS+1);
results_region_based_VD{1,1}    = 'subject';
results_region_based_error{1,1} = 'subject';
for i = 1:nr_VOIS
    results_region_based_VD{1,i+1}    = char(VOIS{i,2});
    results_region_based_error{1,i+1} = char(VOIS{i,2});
end 

% determine the path of the prior data of SPM
tmp             = which('spm.m');
[spm_pth,~,~]   = fileparts(tmp);
spm_pth_priors  = fullfile(spm_pth,'tpm');
% define the brain mask
brain_mask_file = fullfile(spm_pth_priors,'mask_ICV.nii');

curdir = pwd;
for subj = 1:nr_subjects
    clear subjectdir subjectname dir_deriv dir_pet dir_deriv_tmp dir_deriv_anat2 
    clear name_logfile fid go go1 go2 go3 go4 go5 go6 go7
   
    subjectname = SUBJECTS{subj};
    results_region_based_VD{subj+1,1}    = subjectname;
    results_region_based_error{subj+1,1} = subjectname;
    subjectdir  = fullfile(maindir,subjectname);
    if ~isempty(sessiondir)
       dir_pet   = fullfile(fullfile(subjectdir,sessiondir),['pet_' infostring_tracer]);
       dir_deriv = fullfile(derivpetdir,subjectname,sessiondir);
       dir_deriv_tmp   = fullfile(dir_deriv,'tmp');
       dir_deriv_anat2 = fullfile(dir_deriv,'anat');
       dir_firstlevel = fullfile(firstlevelpetdir,subjectname,sessiondir);
        if ~exist(dir_firstlevel,'dir')
            mkdir(dir_firstlevel);
        end
    else
       dir_pet  = fullfile(subjectdir,['pet_' infostring_tracer]); 
       dir_deriv = fullfile(derivpetdir,subjectname);
       dir_deriv_tmp   = fullfile(dir_deriv,'tmp');
       dir_deriv_anat2 = fullfile(dir_deriv,'anat');
       dir_firstlevel = fullfile(firstlevelpetdir,subjectname);
        if ~exist(dir_firstlevel,'dir')
            mkdir(dir_firstlevel);
        end
    end

    
    fprintf('working on subject %s \n',subjectname);
    name_logfile = fullfile(dir_deriv,'LCN12_PET_TSPO_DPA714_log.txt');
    fid  = fopen(name_logfile,'a+');
    fprintf(fid,'subject = %s\n',subjectname);
    fprintf(fid,'LCN12_PET_TSPO_DPA714.m \n');
    fprintf(fid,'%c','-'*ones(1,30));
    fprintf(fid,'\n');
    fprintf(fid,'Processing started at %s\n',datetime('now'));
    fprintf(fid,'\n');
    fprintf(fid,'Settings\n');
    fprintf(fid,'brain_mask_file = %s\n',brain_mask_file);
    fprintf(fid,'logan_start_time = %i min \n',logan_start_time);
    fprintf(fid,'additional_smooth_parametric = %i mm \n',additional_smooth_parametric);
    
    % check if we can find all data in MNI space
    %-------------------------------------------
    [filename_GM,go1]   = LCN_check_filename(dir_deriv_anat2,['wc1' subjectname '*.nii']);
    [filename_WM,go2]   = LCN_check_filename(dir_deriv_anat2,['wc2' subjectname '*.nii']);
    [filename_CSF,go3]  = LCN_check_filename(dir_deriv_anat2,['wc3' subjectname '*.nii']);
    [filename_wPET,go4] = LCN_check_filename(dir_deriv,['w' subjectname '*_' infostring_tracer '*_' infostring_PET '.nii']); 
    [frame_definition_file,go5] = LCN_check_filename(dir_pet,[subjectname '_*' infostring_tracer '_' infostring_PET '_' infostring_frames '.m']); 
    % if we require input function
    if go_input == 1
       [input_file,go6]    = LCN_check_filename(dir_pet,[subjectname '_*' infostring_tracer '_' infostring_input '.m']); 
       [metab_file,go7]    = LCN_check_filename(dir_pet,[subjectname '_*' infostring_tracer '_' infostring_metab '.m']); 
    end
    go = go1.*go2.*go3.*go4.*go5;
    if go_input == 1
       go = go.*go6.*go7;
    end
    
    if go == 1 % all data are found   
       % read the frame definition file
       %+++++++++++++++++++++++++++++++
       fprintf(fid,'Timing file:  %s\n',frame_definition_file);
       copyfile(frame_definition_file,'tmp_frames.m');
       eval('tmp_frames'); % the variable frames_timing is now known
       delete('tmp_frames.m');
       nr_frames = size(frames_timing,1);
       % calculate midscan times (in minutes)
       acqtimes       = frames_timing(:,1:2)/60; % in minutes
       FRAMEDURATION  = acqtimes(:,2)-acqtimes(:,1);
       MIDSCANTIMES   = acqtimes(:,1) + FRAMEDURATION/2;
       
       % reading data in MNI space
       %++++++++++++++++++++++++++
       fprintf('reading data of subject %s\n',subjectname);
       fprintf(fid,'\n');
       fprintf(fid,'Reading data in MNI space started: %s \n',datetime('now'));
       Vref = spm_vol(filename_wPET); 
       % read dynamic data
       dydata = zeros(Vref(1).dim(1),Vref(1).dim(2),Vref(1).dim(3),nr_frames);
       for frame = 1:nr_frames
           tmp = LCN12_read_image([filename_wPET ',' num2str(frame)],Vref(1));
           dydata(:,:,:,frame) = tmp/1000; %in kBq/ml
       end
       fprintf('\n');
       % read segmentations
       fprintf(fid,'GM:  %s\n',filename_GM);
       GMimg = LCN12_read_image(filename_GM,Vref(1));
       fprintf(fid,'WM:  %s\n',filename_WM);
       WMimg = LCN12_read_image(filename_WM,Vref(1));
       fprintf(fid,'CSF: %s\n',filename_CSF);
       CSFimg = LCN12_read_image(filename_CSF,Vref(1));
       fprintf('\n');
       
       if region_based_analysis == 1
          for i = 1:nr_VOIS
              clear filename voi_img
              filename = char(VOIS{i,1});
              VOIname  = char(VOIS{i,2});
              fprintf(fid,'VOI %i: %s\n',i,filename);
              voi_img = LCN12_read_image(filename,Vref(1));
              eval(['VOI_' num2str(i) ' = voi_img;']);
              eval(['VOIname_' num2str(i) ' = VOIname;']);
          end
       end

       % read brain mask in MNI
       fprintf(fid,'brain mask: %s\n',brain_mask_file);
       brain_mask_img = LCN12_read_image(brain_mask_file,Vref(1));
    
       fprintf('reading data of subject %s ... done\n',subjectname);
       fprintf(fid,'reading data ended at %s\n',datetime('now'));
       
       tmp = spm_imatrix(Vref(1).mat);   % tmp(7:9) are the voxel sizes
       voxelsize_wPET = abs(tmp(7:9));
       
       % read the input data and metab data
       %+++++++++++++++++++++++++++++++++++
       fprintf(fid,'Input file:  %s\n',input_file);
       copyfile(input_file,'tmp_input_file.m');
       eval('tmp_input_file'); % the variable input and calibrationfactor_wellcounter are now known
       delete('tmp_input_file.m');
       fprintf(fid,'Metab file:  %s\n',metab_file);
       copyfile(metab_file,'tmp_metab_file.m');
       eval('tmp_metab_file'); % the variable data_metab is now known
       delete('tmp_metab_file.m');

       % analyse the metab data
       %-----------------------
       TIME_METAB = data_metab(:,1)/60; % time p.i. in min
       FRACTION_INTACT_TRACER = data_metab(:,2)/100; % fraction intact tracer
       if size(data_metab,2) == 3
          WEIGHTS_METAB = data_metab(:,3);
       else
          WEIGHTS_METAB = ones(size(data_metab,1),1);
       end
       % normalize WEIGHTS
       WEIGHTS_METAB = WEIGHTS_METAB/sum(WEIGHTS_METAB);        
       % apply a hill fit for the metabolites
       nr_params_hill = length(p0_hill);      
       fine_time  = 0:0.01:max(TIME_METAB);
       nr_samples = length(TIME_METAB);       
       clear pfit res intact_fine
       % fit hill model
       pfit  = fminsearch('LCN_cost_intact_tracer_hill',p0_hill);   
       % determine the error of the fit
       [res] = LCN_cost_intact_tracer_hill(pfit);   
       intact_fine = LCN_calc_intact_tracer_hill(pfit,fine_time);     
       if figures_on == 1 || save_figures == 1
          close('all');
          hfig1 = figure(1);
          set(hfig1,'Name',subjectname);
          plot(TIME_METAB,FRACTION_INTACT_TRACER,'o')
          axis([0 1.1*max(TIME_METAB) 0 1])
          hold on
          plot(fine_time,intact_fine)
          title('hill');
          xlabel('time (min)')
          ylabel('fraction intact tracer')
          if save_figures == 1
             saveas(hfig1,['fig_metab_' subjectname '.fig']);
          else 
             pause
          end
       end

       % analyse the input data (plasma and blood)
       %------------------------------------------
       % calculate corrected inputcurve
       time_input      = input(:,1)/60; % time in min
       [correction]    = LCN_calc_intact_tracer_hill(pfit,time_input);
       uncor_input     = calibrationfactor_wellcounter.*input(:,2)/1000; % in kBq/ml
       corrected_input = calibrationfactor_wellcounter.*correction.*input(:,2)/1000; % in kBq/ml
       uncor_blood     = calibrationfactor_wellcounter.*input(:,3)/1000; % in kBq/ml
       if figures_on == 1 || save_figures == 1
          hfig2 = figure(2);
          set(hfig2,'Name',subjectname);
          subplot(1,2,1);
          plot(time_input,uncor_input,'o');
          hold on
          plot(time_input,corrected_input);
          axis([0 1.1*max(time_input) 0 max(uncor_input)]);
          title('plasma input function');
          xlabel('time (min)');
          ylabel('kBq/ml');
          subplot(1,2,2);
          plot(time_input,uncor_blood,'*r');
          axis([0 1.1*max(time_input) 0 max(uncor_blood)]);
          title('whole blood');
          xlabel('time (min)');
          ylabel('kBq/ml');
          if save_figures == 1
             saveas(hfig2,['fig_input_' subjectname '.fig']);
          else 
             pause
          end
       end
       % determine masks
       brain_mask = (brain_mask_img > 0.5);
       mask = brain_mask;
       
       % load the file excluded_frames
       [filename_excluded_frames,~] = LCN_check_filename(dir_deriv,'excluded_frames_*');
       load(filename_excluded_frames); % now the variable excluded_frames is know
       frames_ok = setdiff([1:nr_frames],excluded_frames);
       
       idx_frames = MIDSCANTIMES > logan_start_time;
       
       if size(frames_ok,2) < (size(FRAMEDURATION,1) - sum(idx_frames) + 2) % we need at least 3 points to fit the Logan regression line
           fprintf('\nWARNING: all frames after logan start time %d are excluded based on head motion, Logan cannot be fit - skipping subject %s\n', logan_start_time, subjectname);
           continue
       end
       
       if region_based_analysis == 1 % perform the region based analysis
          hfig3 = figure(3);
          hfig4 = figure(4);
          for i = 1:nr_VOIS
              clear VD error VOIimg VOIname outputname_VOI
              eval(['VOIimg = VOI_' num2str(i) ';']);
              eval(['VOIname = VOIname_' num2str(i) ';']);
              for frame = 1:nr_frames
                  clear tmp values
                  tmp = dydata(:,:,:,frame);
                  values = tmp((VOIimg > 0.5).*mask > 0);
                  C_MEASURED(frame) = nanmean(values);
              end
              close(3);
              close(4);
              hfig3 = figure(3);
              plot(MIDSCANTIMES,C_MEASURED,'ob:');
              hold on
              plot(MIDSCANTIMES(excluded_frames),C_MEASURED(excluded_frames),'*r:'); % indicating unreliable data due to head movements
              xlabel('time in min');
              ylabel('kBq/ml');
              title(['TAC ' VOIname]);
              hfig4 = figure(4);
              [VD,error] = LCN_LOGAN(MIDSCANTIMES(frames_ok),C_MEASURED(frames_ok),time_input,corrected_input,logan_start_time,4);
              if save_figures == 1
                 saveas(hfig3,['fig_VOI_' VOIname '.fig']);
                 saveas(hfig4,['fig_Logan_VOI_' VOIname '.fig']);
              else 
                 pause
              end
              results_region_based_VD{1+subj,i+1}    = VD;
              results_region_based_error{1+subj,i+1} = error;
          end
       end
       if voxel_based_analysis == 1 % perform the voxel-based analysis
          fprintf('start additional smoothing at %s \n',datetime('now'));
          % smooth the images
          %------------------
          sdydata = zeros(size(dydata));
          for i = 1:nr_frames
              clear tmp stmp
              tmp = squeeze(dydata(:,:,:,i));
              stmp = zeros(size(tmp));
              spm_smooth(tmp,stmp,[additional_smooth_parametric additional_smooth_parametric additional_smooth_parametric]./voxelsize_wPET,0);
              sdydata(:,:,:,i) = stmp;
          end
          fprintf('additional smoothing ended at %s \n',datetime('now'));
          DV_Logan       = zeros(Vref(1).dim(1),Vref(1).dim(2),Vref(1).dim(3));
          DV_Logan_error = zeros(Vref(1).dim(1),Vref(1).dim(2),Vref(1).dim(3));
          fprintf('voxel based analysis started at %s \n',datetime('now'));
          for x = 1:Vref(1).dim(1)
              for y = 1:Vref(1).dim(2)
                  for z = 1:Vref(1).dim(3)
                      if brain_mask(x,y,z) == 1
                         C_MEASURED    = squeeze(sdydata(x,y,z,:));
                         [VD,error] = LCN_LOGAN(MIDSCANTIMES(frames_ok),C_MEASURED(frames_ok),time_input,corrected_input,logan_start_time,0); 
                         DV_Logan(x,y,z)       = VD;
                         DV_Logan_error(x,y,z) = error;
                      end
                  end
              end    
          end
          fprintf('voxel based analysis ended at %s \n',datetime('now'));
          outputname2 = fullfile(dir_firstlevel,['w' subjectname '_DV_Logan.nii']);
          outputname3 = fullfile(dir_firstlevel,['w' subjectname '_DV_Logan_error.nii']);
          Vout1 = LCN12_write_image(DV_Logan,outputname2,'DV Logan',Vref(1).dt(1),Vref(1));        
          Vout2 = LCN12_write_image(DV_Logan_error,outputname3,'DV Logan error',Vref(1).dt(1),Vref(1));        
       end
    else
       fprintf('not all data are found for subject %s. Skipping this subject \n',subjectname); 
       fprintf(fid,'not all data are found for subject %s. Skipping this subject \n',subjectname); 
    end
    fprintf(fid,'++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n\n');
    fclose(fid);
end

% save the results for the region based analysis
if region_based_analysis == 1
    
    % LVO: xlswrite function is deprecated in Matlab and threw an error 
%    xlswrite(outputfile_excel_region_based,results_region_based_VD,'VD');
%    xlswrite(outputfile_excel_region_based,results_region_based_error,'error');   

    fortable_VD = results_region_based_VD(2:end,:);
    table_results_region_based_VD = cell2table(fortable_VD);
    table_results_region_based_VD.Properties.VariableNames = results_region_based_VD(1,:);
    writetable(table_results_region_based_VD, outputfile_excel_region_based, 'Sheet', 'VD');
    
    fortable_error = results_region_based_error(2:end,:);
    table_results_region_based_error = cell2table(fortable_error);
    table_results_region_based_error.Properties.VariableNames = results_region_based_error(1,:);
    writetable(table_results_region_based_error, outputfile_excel_region_based, 'Sheet', 'error');
end