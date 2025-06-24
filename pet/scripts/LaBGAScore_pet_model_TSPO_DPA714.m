% LaBGAScore_pet_model_TSPO_DPA714.m
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
% author: Patrick Dupont, Lukas Van Oudenhove, Lixin Qiu
% date: October 2023
% history: February 2024: an excel file is written which contains for all 
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
%          November 2024: removal of Global variables and adding parallel processing
%          January 2025:  weights for excluded frames set to zero and decay
%                         between measurement samples and start scan is now 
%                         taken into account
%          April 2025:    Logan is also done voxel based
%          May 2025:      Built in option to work with canlab atlas objects (LVO)
%                         Final adaptations to fit LaBGAS file organisation (LVO)
%                         Added atlas labels automatically to output files (LVO)
%          June 2025:     Built in options to choose between voxel-, parcel-wise and ROI analysis (LVO)
%                         Added atlas labels automatically to figure titles/names (LVO)
%                         REMAINS TO BE TESTED
%
% THIS IS RESEARCH SOFTWARE
%__________________________________________________________________________
% @(#)LCN12_PET_TSPO_DPA714.m       v1.4          last modified: 2025/04/30

clear
close all

do_Logan_voxelwise = true;
analysis_type = 'parcel';               % 'parcel' or 'roi'

discoverie_prep_s0_define_directories;

maindir           = BIDSdir;            % directory where the folders of each subject can  be found
sessiondir        = '';                 % if empty, we assume that there is no folder session and the folders anat and pet are directly under the subject folder
infostring_tracer = 'trc-DPA714';
infostring_PET    = 'rec-acdyn_pet';    % the PET data are "subjectname"_"sessiondir"_"infostring_tracer"_"infostring_PET".nii (example =  sub-test_trc-DPA714_rec-acdyn_pet.nii)
infostring_frames = 'frames';           % in the folder pet, we assume a .m file name "subjectname"_"sessiondir"_"infostring_tracer"_"infostring_PET"_"infostring_frames".m (example sub-test_trc-DPA714_rec-acdyn_pet_frames.m)
% infostring_input and infostring_metab need only be defined if we use
% models which require an arterial input function
infostring_input = 'data_blood';        % in the folder pet, we assume a .m file name "subjectname"_"sessiondir"_"infostring_tracer"_"infostring_input".m (example sub-test_trc-DPA714_rec-acdyn_data_blood.nii)
infostring_metab = 'data_metab';        % in the folder pet, we assume a .m file name "subjectname"_"sessiondir"_"infostring_tracer"_"infostring_metab".m (example sub-test_trc-DPA714_rec-acdyn_data_metab.nii)

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

% dir_derivpet = dir(fullfile(derivpetdir,'sub-*'));
dir_derivpet = dir(fullfile(derivpetdir,'sub-*'));
SUBJECTS = {dir_derivpet(:).name}';

% Patrick's code to enter subjects manually
% USE THIS FOR ANALYZING SELECTED SUBJECTS

% SUBJECTS = {
% 'sub-KUL085'
% 'sub-KUL087'
% };

figures_on   = 0;           % if 1, we show the figures for each region and pause. You need to hit a key to proceed.
save_figures = 1;           % if 1, we save the figures to file and pausing of figures is overruled.

% Lukas' code to flexibly load canlab atlas objects for parcel- and
% roi-based analysis

switch analysis_type
    
    case 'parcel'   %  option a - atlas name from load_atlas.m for whole-brain parcellation
        
        atlas_name = 'canlab2024_fine_2mm';
        atlas = load_atlas(atlas_name);
        atlas = downsample_parcellation(atlas,'labels_3'); % downsample to intermediate granularity level, 246 parcels
        atlas.probability_maps = [];
        atlas_labels = atlas.labels;
        atlas_filename = fullfile(secondlevelpetmaskdir,[atlas_name '.nii']);
        if ~isfile(atlas_filename)
            write(atlas,'fname',atlas_filename);
        end
        
    case 'roi'      %  option b - combined roi .nii file generated by https://github.com/labgas/LaBGAScore/blob/main/atlas_mask_tools/LaBGAScore_atlas_rois_from_atlas.m for roi analysis
        
        atlas_name = 'combined_TSPO';
        atlas_name_mat = 'pet_trc-DPA714_combinedTSPO'; % lixin May26 2025
        atlas_mat = fullfile(secondlevelpetmaskdir,[atlas_name_mat '.mat']); % lixin May26 2025
        atlas = load(atlas_mat);
        atlas_labels = atlas.labels;
        atlas_filename = fullfile(secondlevelpetmaskdir,[atlas_name '.nii']);
        
end


% specify the output (we will extend the name with the name of the model
% and write is as an excel file)
outputfile_excel_region_based = fullfile(secondlevelpetresultsdir,['results_' infostring_tracer '_' atlas_name '_VOIS.xlsx']);

% if you want to intersect the VOI region with a subject specific GM or WM map, you need to specify the variables below 
% the same will be done for the parametric image with Logan_input in which case the voxels withing GM or WM after thresholding these maps will be calculated.
intersect_GM = 0;
intersect_WM = 0;
GM_CUTOFF = 0.3;
WM_CUTOFF = 0.3;

additional_smooth_parametric = 6; % kernel size in mm; isotropic Gaussian 3D smoothing; additional smoothing for a voxel based analysis. If GM of WM intersection is selected, smoothing is within this mask

% nr_parpool = 12; % LaBGAS server default
nr_parpool = 20; 

% we will use a multigrid search for the optimal parameters for the rate
% constants for the 2T4k_vb model. Keep in mind that you have
% nr_starting_values_per_dimension^4 number of starting values
% 3 => 81 starting values
% 4 => 256 starting values
% 5 => 625 starting values
% the time for the fitting is roughly proportional to the number of
% starting values.
% the starting value for Vb is fixed but it will be fitted.
Vb0 = 0.05;
nr_starting_values_per_dimension = 3;

%++++ END OF SETTINGS - DO NOT CHANGE BELOW THIS LINE +++++++++++++++++++++

logan_start_time   = 31; % in min (ref Van Weehaeghe et al. J Nucl Med. 2020 Apr;61(4):604-607)
p0_hill        = [50 -1];      % see LCN_calc_intact_tracer_hill for details
k0_2T4k_Vb_lower_bound      = [0.001 0.001 0.001 0.001 0.001];
k0_2T4k_Vb_upper_bound      = [1 1 1 1 1]; 
options = optimoptions('fmincon','Display','off');
STEP = 0.01; % step size for the calculation of integrals in min
% CALC_OPTION: parameter determing the way we calculate the output of a
% model.
%   1 - the model is calculated as the integral of the output concentration
%       devided by the frameduration.
%   2 - the model is calculated at the midscantime.
CALC_OPTION = 2;
window_size = 3; % for temporal median filtering
min_number_voxels = 5; % minimum number of voxels in the parcel that are used for the calculation of the TAC.
thalf_F18 = 1.82871*60; % half life in min of 18F (REF: García-Toraño E, Medina VP, Ibarra MR. The half-life of 18F. Appl Radiat Isot. (2010) 68(7-8):1561-5)
%--------------------------------------------------------------------------
curdir = pwd;
% determine the path of the prior data of SPM
tmp             = which('spm.m');
[spm_pth,~,~]   = fileparts(tmp);
spm_pth_priors  = fullfile(spm_pth,'tpm');
% define the brain mask
brain_mask_file = fullfile(spm_pth_priors,'mask_ICV.nii');

nr_subjects = size(SUBJECTS,1);
% nr_models   = size(model_list,1);
%% ATLAS PARCELS
if exist('atlas_filename','var') == 1 && ~isempty(atlas_filename) == 1
   % read atlas
   [atlas_orig,Vatlas_orig] = LCN12_read_image(atlas_filename);
   atlas_orig = round(atlas_orig);
   parcel_values = setdiff(unique(atlas_orig(:)),0); % assuming 0 is background
   % check if there is a .m file with the same name as the atlas
   % with the VOIdetails such as the name
   [pth,name,ext] = fileparts(atlas_filename);
   if exist(fullfile(pth,[name '.m']),'file') == 2
      % execute this m-file
      curdir = pwd;
      cd(pth);
      eval(name); % now VOIdetails is known, a n x 2 cell array with each row the value of a parcel in the atlas and the corresponding parcel name
      cd(curdir);
      parcel_values_atlas = cell2mat(VOIdetails(:,1));
      parcel_names = VOIdetails(:,2);
   else
      parcel_values_atlas = parcel_values;
      parcel_names = char(atlas_labels);
   end
   nr_VOIS = length(parcel_values);
   % initialize output
   varNames = cell(1,nr_VOIS+1);
   for voi = 1:nr_VOIS
       index_value = find(parcel_values_atlas == parcel_values(voi));
       varNames(1,voi+1) = {[num2str(parcel_values(voi)) '-' char(parcel_names(index_value,:))]};
   end
else
   fprintf('You need to define an atlas \n');
   return;
end
% initialize output
sz = [nr_subjects nr_VOIS+1];
varNames{1}  = 'subject';
varTypes(1,1)           = {'string'};
varTypes(1,2:nr_VOIS+1) = {'double'};
results_region_based_Logan_input_DV    = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames);
results_region_based_Logan_input_error = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames);
results_region_based_2T4k_Vb_K1    = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames);
results_region_based_2T4k_Vb_k2    = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames);
results_region_based_2T4k_Vb_k3    = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames);
results_region_based_2T4k_Vb_k4    = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames);
results_region_based_2T4k_Vb_Vb    = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames);
results_region_based_2T4k_Vb_DV    = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames);
results_region_based_2T4k_Vb_error = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames);

% determing the grid of starting values within the search space
%--------------------------------------------------------------
k0_all_values = zeros(nr_starting_values_per_dimension^4,5);
counter = 0;
for d1 = 1:nr_starting_values_per_dimension
    val1 = k0_2T4k_Vb_lower_bound(1) + d1*(k0_2T4k_Vb_upper_bound(1) - k0_2T4k_Vb_lower_bound(1))/(nr_starting_values_per_dimension+1);
    for d2 = 1:nr_starting_values_per_dimension
        val2 = k0_2T4k_Vb_lower_bound(2) + d2*(k0_2T4k_Vb_upper_bound(2) - k0_2T4k_Vb_lower_bound(2))/(nr_starting_values_per_dimension+1);
        for d3 = 1:nr_starting_values_per_dimension
            val3 = k0_2T4k_Vb_lower_bound(3) + d3*(k0_2T4k_Vb_upper_bound(3) - k0_2T4k_Vb_lower_bound(3))/(nr_starting_values_per_dimension+1);
            for d4 = 1:nr_starting_values_per_dimension
                val4 = k0_2T4k_Vb_lower_bound(4) + d4*(k0_2T4k_Vb_upper_bound(4) - k0_2T4k_Vb_lower_bound(4))/(nr_starting_values_per_dimension+1);
                for d5 = 1:nr_starting_values_per_dimension
                    val5 = Vb0;
                    counter = counter + 1;
                    k0_all_values(counter,:) = [val1 val2 val3 val4 val5]; 
                end
            end
        end
    end
end
nr_k0 = size(k0_all_values,1);

%% loop over the subjects
parpool(nr_parpool);
for subj = 1:nr_subjects
    
    clear subjectdir subjectname dir_anat dir_pet dir_deriv* dir_firstlevel dir_tmp dir_anat2 DV_Logan
    clear tmp_frames frames_timing nr_frames acqtimes FRAMEDURATION MIDSCANTIMES WEIGHTS_PET
    clear name_logfile fid go go_GM go_WM go_CSF go_PET go_frames go_input go_metab Vref Vref_tmp
    
    subjectname = SUBJECTS{subj};
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

    if save_figures == 1
       dir_figures = fullfile(dir_deriv,'figures');
       if exist(dir_figures,'dir') ~= 7
          mkdir(dir_figures);
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
    fprintf(fid,'figures_on    = %i\n',figures_on);
    fprintf(fid,'save_figures  = %i\n',save_figures);
    fprintf(fid,'atlas_filename = %s\n',atlas_filename);
    if exist('intersect_GM','var') == 1
       fprintf(fid,'intersect_GM = %i\n',intersect_GM);
    else
       fprintf(fid,'intersect_GM = not defined \n');
    end
    if exist('intersect_WM','var') == 1
       fprintf(fid,'intersect_WM = %i\n',intersect_WM);
    else
       fprintf(fid,'intersect_WM = not defined \n');
    end
    if exist('GM_CUTOFF','var') == 1
       fprintf(fid,'GM_CUTOFF = %4,2f\n',GM_CUTOFF);
    else
       fprintf(fid,'GM_CUTOFF = not defined \n');
    end
    if exist('WM_CUTOFF','var') == 1
       fprintf(fid,'WM_CUTOFF = %4.2f\n',WM_CUTOFF);
    else
       fprintf(fid,'WM_CUTOFF = not defined \n');
    end
    fprintf(fid,'logan_start_time (min) = %i \n',logan_start_time);
    fprintf(fid,'p0_hill = [%4.2f %4.2f] \n',p0_hill(1),p0_hill(2));
    fprintf(fid,'k0_2T4k_Vb_lower_bound      = [%4.2f %4.2f %4.2f %4.2f %4.2f] \n',k0_2T4k_Vb_lower_bound(1),k0_2T4k_Vb_lower_bound(2),k0_2T4k_Vb_lower_bound(3),k0_2T4k_Vb_lower_bound(4),k0_2T4k_Vb_lower_bound(5));
    fprintf(fid,'k0_2T4k_Vb_upper_bound      = [%4.2f %4.2f %4.2f %4.2f %4.2f] \n',k0_2T4k_Vb_upper_bound(1),k0_2T4k_Vb_upper_bound(2),k0_2T4k_Vb_upper_bound(3),k0_2T4k_Vb_upper_bound(4),k0_2T4k_Vb_upper_bound(5));
        
    % check if we can find all data in MNI space
    [filename_GM,go_GM]    = LCN_check_filename(dir_deriv_anat2,['wc1' subjectname '*.nii']);
    [filename_WM,go_WM]    = LCN_check_filename(dir_deriv_anat2,['wc2' subjectname '*.nii']);
    [filename_CSF,go_CSF]  = LCN_check_filename(dir_deriv_anat2,['wc3' subjectname '*.nii']);
    [filename_wPET,go_PET] = LCN_check_filename(dir_deriv,['w' subjectname '*_' infostring_tracer '*_' infostring_PET '.nii']); 
    [frame_definition_file,go_frames] = LCN_check_filename(dir_pet,[subjectname '_*' infostring_tracer '_' infostring_PET '_' infostring_frames '.m']); 
    [input_file,go_input]  = LCN_check_filename(dir_pet,[subjectname '_*' infostring_tracer '_' infostring_input '.m']); 
    [metab_file,go_metab]  = LCN_check_filename(dir_pet,[subjectname '_*' infostring_tracer '_' infostring_metab '.m']); 

    if go_frames == 1 % frame defintion file found   
       % read the frame definition file
       fprintf(fid,'Timing file:  %s\n',frame_definition_file);
       copyfile(frame_definition_file,'tmp_frames.m');
       tmp_frames; % the variable frames_timing is now known
       delete('tmp_frames.m');
       nr_frames = size(frames_timing,1);
       % calculate midscan times (in minutes)
       acqtimes       = frames_timing(:,1:2)/60; % in minutes
       FRAMEDURATION  = acqtimes(:,2)-acqtimes(:,1);
       MIDSCANTIMES   = acqtimes(:,1) + FRAMEDURATION/2;
    else
       fprintf(fid,'frame definition file not found for subject %s \n',subjectname);
       continue;
    end
    if go_input*go_metab == 1
       clear ref_time scanner_time
       % read the input data and metab data
       fprintf(fid,'Input file:  %s\n',input_file);
       copyfile(input_file,'tmp_input_file.m');
       tmp_input_file; % the variable input and calibrationfactor_wellcounter are now known
       delete('tmp_input_file.m');
       fprintf(fid,'Metab file:  %s\n',metab_file);
       copyfile(metab_file,'tmp_metab_file.m');
       tmp_metab_file; % the variable data_metab is now known
       delete('tmp_metab_file.m');
       % analyse the metab data
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
       % nr_params_hill = length(p0_hill);      
       % nr_samples = length(TIME_METAB);       
       clear pfit res intact_fine
       % fit hill model
       cost_function = @(p) norm(WEIGHTS_METAB .* (LCN_calc_intact_tracer_hill(p,TIME_METAB)' - FRACTION_INTACT_TRACER));
       [pfit, error] = fminsearch(cost_function, p0_hill);       
       if figures_on == 1 || save_figures == 1
          close('all');
          hfig1 = figure(1);
          set(hfig1,'Name',subjectname);
          plot(TIME_METAB,FRACTION_INTACT_TRACER,'o')
          axis([0 1.1*max(TIME_METAB) 0 1])
          hold on
          fine_time   = (0:0.01:max(TIME_METAB))';
          intact_fine = LCN_calc_intact_tracer_hill(pfit,fine_time);     
          plot(fine_time,intact_fine)
          title('hill');
          xlabel('time (min)')
          ylabel('fraction intact tracer')
          if save_figures == 1
             cd(dir_figures);
             saveas(hfig1,['fig_metab_' subjectname '.fig']);
          else 
             pause
          end
       end
       % analyse the input data (plasma and blood)
       % if the variables ref_time and scanner_time exist, we need to
       % correct for decay between measurement of samples and start scan.
       % if they are not available, it means that the data are not yet
       % corrected for this time difference and it will be done here.
       if exist('ref_time','var') == 1 && exist('scanner_time','var') == 1
          % correct blood data for decay between start scan and start measurement and well counter
          dt = minutes(duration(datetime(ref_time)-datetime(scanner_time),'Format','s'));
          decay_factor = exp(log(2)*dt/thalf_F18);
          input(:,2) = input(:,2).*decay_factor;
          input(:,3) = input(:,3).*decay_factor;
       end

       % calculate corrected inputcurve
       time_input      = input(:,1)/60; % time in min
       [correction]    = LCN_calc_intact_tracer_hill(pfit,time_input);
       uncor_input     = calibrationfactor_wellcounter.*input(:,2)/(60*1000); % in kBq/ml
       corrected_input = calibrationfactor_wellcounter.*correction.*input(:,2)/(60*1000); % in kBq/ml
       uncor_blood     = calibrationfactor_wellcounter.*input(:,3)/(60*1000); % in kBq/ml
       TIME            = [0; time_input];
       FINETIMES       = (0:STEP:max(TIME))';
       CA              = [0; corrected_input];
       CA_WB           = [0; uncor_blood];
       CA_FINETIMES    = interp1q(TIME,CA,FINETIMES);
       CA_MIDSCANTIMES = interp1q(TIME,CA,MIDSCANTIMES);
       CA_WB_FINETIMES = interp1q(TIME,CA_WB,FINETIMES);
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
             cd(dir_figures);
             saveas(hfig2,['fig_input_' subjectname '.fig']);
          else 
             pause
          end
       end
    end
    
    if go_PET == 1 % dynamic data are found   
       % reading data in MNI space
       fprintf('reading data of subject %s',subjectname);
       fprintf(fid,'\n');
       fprintf(fid,'Reading data in MNI space started: %s \n',datetime('now'));
       Vref_tmp = spm_vol(filename_wPET); 
       Vref = Vref_tmp(1);
       % read dynamic data
       dydata = zeros(Vref.dim(1),Vref.dim(2),Vref.dim(3),nr_frames);
       for frame = 1:nr_frames
           clear tmp
           tmp = LCN12_read_image([filename_wPET ',' num2str(frame)],Vref);
           dydata(:,:,:,frame) = tmp/1000; %in kBq/ml
       end
       % load the file excluded_frames
       [filename_excluded_frames,~] = LCN_check_filename(dir_deriv,'excluded_frames_*pet.mat');
       load(filename_excluded_frames); % now the variable excluded_frames is know
       frames_ok = setdiff(1:nr_frames,excluded_frames)';
       
       idx_frames = MIDSCANTIMES > logan_start_time;
       
       if size(frames_ok,1) < (size(FRAMEDURATION,1) - sum(idx_frames) + 2) % we need at least 3 points to fit the Logan regression line
           fprintf('\nWARNING: all frames after logan start time %d are excluded based on head motion, Logan cannot be fit - skipping subject %s\n', logan_start_time, subjectname);
           continue
       end
       
       % read brain mask in MNI
       fprintf(fid,'brain mask: %s\n',brain_mask_file);
       brain_mask_img = LCN12_read_image(brain_mask_file,Vref);
       % determine masks
       brain_mask = (brain_mask_img > 0.5);
       % get the voxel size
       tmp = spm_imatrix(Vref.mat);   % tmp(7:9) are the voxel sizes
       voxelsize_wPET = abs(tmp(7:9));
    else
       fprintf('Dynamic data of subject %s not found \n',subjectname);
       continue;
    end
    if go_GM*go_WM*go_CSF == 1
       % read segmentations
       fprintf(fid,'GM:  %s\n',filename_GM);
       GMimg = LCN12_read_image(filename_GM,Vref);
       fprintf(fid,'WM:  %s\n',filename_WM);
       WMimg = LCN12_read_image(filename_WM,Vref);
       fprintf(fid,'CSF: %s\n',filename_CSF);
       CSFimg = LCN12_read_image(filename_CSF,Vref);
    end
    fprintf('... done\n');
    fprintf(fid,'reading data ended at %s\n',datetime('now'));
    if figures_on == 1 || save_figures == 1
       hfig3 = figure(3);
       hfig4 = figure(4);
    end
    results_region_based_Logan_input_DV(subj,1)    = {subjectname};
    results_region_based_Logan_input_error(subj,1) = {subjectname};
    atlas_Logan_input_DV    = zeros(size(atlas_orig));
    atlas_Logan_input_error = zeros(size(atlas_orig));
    results_region_based_2T4k_Vb_K1(subj,1)    = {subjectname};
    results_region_based_2T4k_Vb_k2(subj,1)    = {subjectname};
    results_region_based_2T4k_Vb_k3(subj,1)    = {subjectname};
    results_region_based_2T4k_Vb_k4(subj,1)    = {subjectname};
    results_region_based_2T4k_Vb_Vb(subj,1)    = {subjectname};
    results_region_based_2T4k_Vb_DV(subj,1)    = {subjectname};
    results_region_based_2T4k_Vb_error(subj,1) = {subjectname};
    atlas_2T4k_Vb_K1    = zeros(size(atlas_orig));
    atlas_2T4k_Vb_k2    = zeros(size(atlas_orig));
    atlas_2T4k_Vb_k3    = zeros(size(atlas_orig));
    atlas_2T4k_Vb_k4    = zeros(size(atlas_orig));
    atlas_2T4k_Vb_Vb    = zeros(size(atlas_orig));
    atlas_2T4k_Vb_VT    = zeros(size(atlas_orig));
    atlas_2T4k_Vb_error = zeros(size(atlas_orig));
    % read the atlas in the matrix of the subject
    atlas_s = LCN12_read_image(atlas_filename,Vref);
    % atlas contains non-integer values in voxels which are a mixture. We focus
    % only on voxels which are of one specific type
    % therefore, set the value of other voxels to 0
    atlas_s(mod(atlas_s,1) > 0) = 0;
    % adapt weight of the frames 
    WEIGHTS_PET = ones(length(MIDSCANTIMES),1);
    WEIGHTS_PET(excluded_frames) = 0;
    WEIGHTS_PET    = WEIGHTS_PET./sum(WEIGHTS_PET);
    fprintf('atlas based analysis started at %s \n',datetime('now'));
    fprintf(fid,'atlas based analysis started at %s \n',datetime('now'));
    
    % read the VOIS
    for voi = 1:nr_VOIS
        clear index_value VOIimg VOIname C_MEASURED
        index_value = find(parcel_values_atlas == parcel_values(voi));
        VOIname  = char(parcel_names(index_value,:));
        fprintf(fid,'VOI %i: %s\n',index_value,VOIname);
        VOIimg = (atlas_s == parcel_values(voi));
        % intersect VOI with subject specific GM or WM
        if intersect_GM == 1
           VOIimg = (VOIimg > 0.5).*(GMimg > GM_CUTOFF);
        elseif intersect_WM  == 1
           VOIimg = (VOIimg > 0.5).*(WMimg > WM_CUTOFF);
        end
        if nnz(VOIimg(:) > 0) < min_number_voxels 
           continue;
        end
        C_MEASURED = zeros(nr_frames,1);
        for frame = 1:nr_frames
            clear tmp values
            tmp = dydata(:,:,:,frame);
            values = tmp((VOIimg > 0.5).*brain_mask > 0);
            C_MEASURED(frame) = mean(values(~isnan(values)));
        end
        % temporal filtering using a median filter
        C_MEASURED = medfilt1(C_MEASURED, window_size);
        if sum(abs(C_MEASURED)) == 0
           continue;
        end
        if sum(isnan(C_MEASURED)) > 0
           continue;
        end
        % calculate all the models
        % LOGAN_INPUT
        clear VD error hfig3
        if figures_on == 1 || save_figures == 1 
           close(3);
           hfig3 = figure(3);
           [DV,error] = LCN_LOGAN(MIDSCANTIMES(frames_ok),C_MEASURED(frames_ok),time_input,corrected_input,logan_start_time,3);
           title(['Logan plot: ' num2str(index_value) '-' LCN_string(VOIname)]);
           % set(hfig4,'NumberTitle','off')
           set(hfig3,'Name',[subjectname '- Logan\_input']);
           if save_figures == 1
              cd(dir_figures);
              saveas(hfig3,['fig_Logan_input_VOI_' num2str(index_value) '-' VOIname '.fig']);
           else 
              pause
           end
        else
           [DV,error] = LCN_LOGAN(MIDSCANTIMES(frames_ok),C_MEASURED(frames_ok),time_input,corrected_input,logan_start_time,0);
        end
        results_region_based_Logan_input_DV(subj,1+voi)    = {DV};
        results_region_based_Logan_input_error(subj,1+voi) = {error};
        atlas_Logan_input_DV(atlas_orig == parcel_values(voi))    = DV;
        atlas_Logan_input_error(atlas_orig == parcel_values(voi)) = error;
        
        % 2T4k_Vb
        clear tmpA tmpB index_min VD error_2T4k_Vb kfit_2T4k_Vb cost_function fit_values VT_new hfig4 axvalues
        % start new fit from different starting positions
        tmpA = zeros(nr_k0,5);
        tmpB = zeros(nr_k0,1);
        cost_function = @(params) norm(WEIGHTS_PET.*(LCN_calc2_model_2T4k_Vb(params,MIDSCANTIMES,FRAMEDURATION,CA_FINETIMES,CA_WB_FINETIMES,FINETIMES,CALC_OPTION,STEP) - C_MEASURED));
        parfor j = 1:nr_k0
               k0 = k0_all_values(j,:);
               [tmpA(j,:),tmpB(j)] = fmincon(cost_function,k0,[],[],[],[],k0_2T4k_Vb_lower_bound,k0_2T4k_Vb_upper_bound,[],options);
        end
        index_min = find(tmpB == min(tmpB));
        kfit_2T4k_Vb = tmpA(index_min(1),:);
        results_region_based_2T4k_Vb_K1(subj,1+voi)    = {kfit_2T4k_Vb(1)};
        results_region_based_2T4k_Vb_k2(subj,1+voi)    = {kfit_2T4k_Vb(2)};
        results_region_based_2T4k_Vb_k3(subj,1+voi)    = {kfit_2T4k_Vb(3)};
        results_region_based_2T4k_Vb_k4(subj,1+voi)    = {kfit_2T4k_Vb(4)};
        results_region_based_2T4k_Vb_Vb(subj,1+voi)    = {kfit_2T4k_Vb(5)};
        results_region_based_2T4k_Vb_DV(subj,1+voi)    = {kfit_2T4k_Vb(1)*(1+kfit_2T4k_Vb(3)/kfit_2T4k_Vb(4))/kfit_2T4k_Vb(2)};
        results_region_based_2T4k_Vb_error(subj,1+voi) = {tmpB(index_min(1))./mean(C_MEASURED)}; % relative error
        atlas_2T4k_Vb_K1(atlas_orig == parcel_values(voi))    = kfit_2T4k_Vb(1);
        atlas_2T4k_Vb_k2(atlas_orig == parcel_values(voi))    = kfit_2T4k_Vb(2);
        atlas_2T4k_Vb_k3(atlas_orig == parcel_values(voi))    = kfit_2T4k_Vb(3);
        atlas_2T4k_Vb_k4(atlas_orig == parcel_values(voi))    = kfit_2T4k_Vb(4);
        atlas_2T4k_Vb_Vb(atlas_orig == parcel_values(voi))    = kfit_2T4k_Vb(5);
        atlas_2T4k_Vb_VT(atlas_orig == parcel_values(voi))    = kfit_2T4k_Vb(1)*(1+kfit_2T4k_Vb(3)/kfit_2T4k_Vb(4))/kfit_2T4k_Vb(2);
        atlas_2T4k_Vb_error(atlas_orig == parcel_values(voi)) = tmpB(index_min(1))./mean(C_MEASURED);
        if figures_on == 1 || save_figures == 1 
           close(4);
           hfig4 = figure(4);
           plot(MIDSCANTIMES,C_MEASURED,'ob:');
           hold on
           plot(MIDSCANTIMES(excluded_frames),C_MEASURED(excluded_frames),'*r:'); % indicating unreliable data due to head movements
           xlabel('time in min');
           ylabel('kBq/ml');
           title(['TAC ' num2str(index_value) '-' LCN_string(VOIname)]);
           for j = 1:nr_k0
               plot(MIDSCANTIMES,LCN_calc2_model_2T4k_Vb(tmpA(j,:),MIDSCANTIMES,FRAMEDURATION,CA_FINETIMES,CA_WB_FINETIMES,FINETIMES,CALC_OPTION,STEP),'b')
           end
           plot(MIDSCANTIMES,LCN_calc2_model_2T4k_Vb(kfit_2T4k_Vb,MIDSCANTIMES,FRAMEDURATION,CA_FINETIMES,CA_WB_FINETIMES,FINETIMES,CALC_OPTION,STEP),'r')
           axvalues = axis;
           text(axvalues(2)*0.5,axvalues(3)+0.1*(axvalues(4)-axvalues(3)),['VD = ' num2str(kfit_2T4k_Vb(1)*(1+kfit_2T4k_Vb(3)/kfit_2T4k_Vb(4))/kfit_2T4k_Vb(2))]);
           % set(hfig4,'NumberTitle','off')
           set(hfig4,'Name',[subjectname '-2T4k\_Vb']);
           if save_figures == 1
              cd(dir_figures);
              saveas(hfig4,['fig_2T4k_Vb_' num2str(index_value) '-' VOIname '.fig']);
           else 
              pause
           end
        end
    end
    fprintf('atlas based analysis ended at %s \n',datetime('now'));
    fprintf(fid,'atlas based analysis ended at %s \n',datetime('now'));

    % write the results to an image in which each parcel has the same value
    %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    % Logan_input
    outputname1 = fullfile(dir_firstlevel,['w' subjectname '_atlas_Logan_input_DV.nii']);
    outputname2 = fullfile(dir_firstlevel,['w' subjectname '_atlas_Logan_input_error.nii']);
    LCN12_write_image(atlas_Logan_input_DV,outputname1,'Logan input DV',64,Vatlas_orig);        
    LCN12_write_image(atlas_Logan_input_error,outputname2,'Logan input error',64,Vatlas_orig);        
    % 2T4k_Vb
    outputname1 = fullfile(dir_firstlevel,['w' subjectname '_atlas_2T4k_Vb_K1.nii']);
    outputname2 = fullfile(dir_firstlevel,['w' subjectname '_atlas_2T4k_Vb_k2.nii']);
    outputname3 = fullfile(dir_firstlevel,['w' subjectname '_atlas_2T4k_Vb_k3.nii']);
    outputname4 = fullfile(dir_firstlevel,['w' subjectname '_atlas_2T4k_Vb_k4.nii']);
    outputname5 = fullfile(dir_firstlevel,['w' subjectname '_atlas_2T4k_Vb_Vb.nii']);
    outputname6 = fullfile(dir_firstlevel,['w' subjectname '_atlas_2T4k_Vb_VT.nii']);
    outputname7 = fullfile(dir_firstlevel,['w' subjectname '_atlas_2T4k_Vb_error.nii']);
    LCN12_write_image(atlas_2T4k_Vb_K1,outputname1,'2T4k_Vb_K1',64,Vatlas_orig);        
    LCN12_write_image(atlas_2T4k_Vb_k2,outputname2,'2T4k_Vb_k2',64,Vatlas_orig);        
    LCN12_write_image(atlas_2T4k_Vb_k3,outputname3,'2T4k_Vb_k3',64,Vatlas_orig);        
    LCN12_write_image(atlas_2T4k_Vb_k4,outputname4,'2T4k_Vb_k4',64,Vatlas_orig);        
    LCN12_write_image(atlas_2T4k_Vb_Vb,outputname5,'2T4k_Vb_Vb',64,Vatlas_orig);        
    LCN12_write_image(atlas_2T4k_Vb_VT,outputname6,'2T4k_Vb_VT',64,Vatlas_orig);        
    LCN12_write_image(atlas_2T4k_Vb_error,outputname7,'2T4k_Vb_error',64,Vatlas_orig);        

    if do_Logan_voxelwise
    
        % generate a voxelbased parametric image for Logan_input
        clear DV_Logan DV_Logan_error outputname1 outputname2 Vout1 Vout2
        dimx = Vref.dim(1);
        dimy = Vref.dim(2);
        dimz = Vref.dim(3);
        fprintf('start additional smoothing at %s \n',datetime('now'));
        fprintf(fid,'start additional smoothing at %s \n',datetime('now'));
        % smooth the images
        sdydata = zeros(dimx*dimy*dimz,size(dydata,4));
        if intersect_GM == 1
           brain_mask = (brain_mask > 0).*(GMimg > GM_CUTOFF);
        end
        if intersect_WM == 1
           brain_mask = (brain_mask > 0).*(WMimg > WM_CUTOFF);
        end
        for frame = 1:nr_frames
            clear tmp stmp
            tmp = squeeze(dydata(:,:,:,frame)).*(brain_mask > 0);
            stmp = LCN12_smooth(tmp,[additional_smooth_parametric additional_smooth_parametric additional_smooth_parametric]./voxelsize_wPET,brain_mask);
            sdydata(:,frame) = stmp(:);
        end
        fprintf('additional smoothing ended at %s \n',datetime('now'));
        fprintf(fid,'additional smoothing ended at %s \n',datetime('now'));
        indexlist = find(brain_mask(:) == 1);
        datavalues = sdydata(indexlist,:);
        tmpA = [];
        tmpB = [];
        tmp1 = zeros(dimx*dimy*dimz,1);
        tmp2 = zeros(dimx*dimy*dimz,1);
        fprintf('voxel based analysis for Logan_input started at %s \n',datetime('now'));
        fprintf(fid,'voxel based analysis for Logan_input started at %s \n',datetime('now'));
        datavalues_ok   = datavalues(:,frames_ok);
        midscantimes_ok = MIDSCANTIMES(frames_ok);
        parfor i = 1:length(indexlist)
               C_MEASURED    = datavalues_ok(i,:)';
               [DV,error] = LCN_LOGAN(midscantimes_ok,C_MEASURED,time_input,corrected_input,logan_start_time,0);
               tmpA(i) = DV;
               tmpB(i) = error;                     
        end
        tmp1(indexlist) = tmpA;
        tmp2(indexlist) = tmpB;
        DV_Logan        = reshape(tmp1,dimx,dimy,dimz);
        DV_Logan_error  = reshape(tmp2,dimx,dimy,dimz);
        outputname1 = fullfile(dir_firstlevel,['w' subjectname '_Logan_input_DV.nii']);
        outputname2 = fullfile(dir_firstlevel,['w' subjectname '_Logan_input_error.nii']);
        Vout1 = LCN12_write_image(DV_Logan,outputname1,'DV Logan',Vref.dt(1),Vref);        
        Vout2 = LCN12_write_image(DV_Logan_error,outputname2,'DV Logan error',Vref.dt(1),Vref);        
        fprintf('voxel based analysis for Logan_input ended at %s \n',datetime('now'));
        fprintf(fid,'voxel based analysis for Logan_input ended at %s \n',datetime('now'));
    
    end
    
end

%% save the results for the region based analysis
%++++++++++++++++++++++++++++++++++++++++++++++++
% Logan_input
outputfile_excel1 = [outputfile_excel_region_based '_Logan_input.xlsx'];
writetable(results_region_based_Logan_input_DV,outputfile_excel1,'Sheet','DV');
writetable(results_region_based_Logan_input_error,outputfile_excel1,'Sheet','error');
% 2T4k_Vb
outputfile_excel2 = [outputfile_excel_region_based '_2T4k_Vb.xlsx'];
writetable(results_region_based_2T4k_Vb_K1,outputfile_excel2,'Sheet','K1');
writetable(results_region_based_2T4k_Vb_k2,outputfile_excel2,'Sheet','k2');
writetable(results_region_based_2T4k_Vb_k3,outputfile_excel2,'Sheet','k3');
writetable(results_region_based_2T4k_Vb_k4,outputfile_excel2,'Sheet','k4');
writetable(results_region_based_2T4k_Vb_Vb,outputfile_excel2,'Sheet','Vb');
writetable(results_region_based_2T4k_Vb_DV,outputfile_excel2,'Sheet','DV');
writetable(results_region_based_2T4k_Vb_error,outputfile_excel2,'Sheet','error');
cd(curdir);
fprintf('ALL DONE\n');