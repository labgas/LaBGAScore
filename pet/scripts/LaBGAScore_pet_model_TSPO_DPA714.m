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
%          Sept 2025:     Integration of LCN12_PET_TSPO_DPA714_compare_models.m (v0.2)
%                         ~ includes 2T4k models with (ir)reversible endothelial binding
%                         component and fit statistics to compare models
%                         (LVO)
%
% THIS IS RESEARCH SOFTWARE
%__________________________________________________________________________
% @(#)LCN12_PET_TSPO_DPA714.m       v2.0          last modified: 2025/09/22

clear
close all

%% PREP WORK, SET INFO AND OPTIONS
%--------------------------------------------------------------------------

% SET DIRECTORIES

LaBGAScore_prep_s0_define_directories; % MAKE STUDY-SPECIFIC

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

    
% SELECT SUBJECTS

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
% 'sub-KUL034'
% };


% SET FIGURE/IMAGE OPTIONS

figures_on   = 0;           % if 1, we show the figures for each region and pause. You need to hit a key to proceed.
save_figures = 1;           % if 1, we save the figures to file and pausing of figures is overruled.
save_parcelimgs = 0;        % if 1, we save parcel-wise images for every parameter in every model for every subject
do_voxelwise_logan = 0;     % if 1, we perform a voxel-wise Logan analysis in addition to the parcel-wise analyses 
pet_space_reference = 0;    % if 1, atlas is read in pet space. Voxels with non-integer values (i.e. voxels at the border of a region) will be excluded. If 0, atlas space is used.
suffix = 'atlas_space';     % use if you want to run multiple models

% DEFINE ATLAS/ROIS

% Lukas' code to flexibly load canlab atlas objects

%  option a - atlas name from load_atlas.m for whole-brain parcellation
atlas_name = 'canlab2024_fine_2mm';
atlas = load_atlas(atlas_name);
atlas = downsample_parcellation(atlas,'labels_3'); % downsample to intermediate granularity level, 246 parcels
atlas.probability_maps = [];
atlas_filename = fullfile(secondlevelpetmaskdir,[atlas_name '.nii']);
if ~isfile(atlas_filename)
    write(atlas,'fname',atlas_filename);
end

%  option b - combined roi .nii file generated by
%  https://github.com/labgas/LaBGAScore/blob/main/atlas_mask_tools/LaBGAScore_atlas_rois_from_atlas.m
%  for roi analysis
% atlas_name = 'combined_inflammation_regions';
% atlas_name = 'combined_TSPO';
% atlas_mat = 'pet_trc-DPA714_combinedTSPO'; % lixin May26 2025
% atlas_filename = fullfile(secondlevelpetmaskdir,[atlas_name '.nii']);
% atlas_mat = fullfile(secondlevelpetmaskdir,[atlas_mat '.mat']); % lixin May26 2025
% load(atlas_mat);


% SPECIFY BASE NAME FOR THE OUTPUT FILES
% we will extend the name with the name of the model and write as an excel file)

outputfile_excel_region_based = fullfile(secondlevelpetresultsdir,['results_' infostring_tracer '_' atlas_name '_VOIS']);


% SET OPTIONS FOR INTERSECTING WITH GRAY MATTER

% if you want to intersect the VOI region with a subject specific GM or WM map, you need to specify the variables below 
% the same will be done for the parametric image with Logan_input in which case the voxels withing GM or WM after thresholding these maps will be calculated.
intersect_GM = 0;
intersect_WM = 0;
GM_CUTOFF = 0.3;
WM_CUTOFF = 0.3;


% SET SMOOTHING KERNEL FOR VOXEL-BASED PARAMETRIC LOGAN IMAGES

additional_smooth_parametric = 6; % kernel size in mm; isotropic Gaussian 3D smoothing; additional smoothing for a voxel based analysis. If GM of WM intersection is selected, smoothing is within this mask


% SET NUMBER OF PARALLEL PROCESSES

nr_parpool = 12; % LaBGAS server default
% nr_parpool = 20; 

%+++++++++++++++ DO NOT CHANGE BELOW THIS LINE ++++++++++++++++++++++++++++


%% SETTINGS
%---------------------- settings ------------------------------------------

min_number_voxels = 25; % minimum number of voxels in the parcel that are used for the calculation of the TAC.
% we will use a multigrid search for the optimal parameters for the rate
% constants for each model. Keep in mind that you have
% nr_starting_values_per_dimension^4 number of starting values for a model
% with 4 rate constants that you vary and
% nr_starting_values_per_dimension^6 for a model with 6 parameters that you
% would like to vary.
% the time for the fitting is roughly proportional to the number of
% starting values.
% the starting value for Vb is fixed but it will be fitted.
Vb0 = 0.05;
nr_starting_values_per_dimension = 3;
model_list = {
    'Logan_input'
    '2T4k'
    '2T4k_Vb'
    '2T4k_vasc1k'
    '2T4k_vasc2k'
};
logan_start_time   = 31; % in min (ref Van Weehaeghe et al. J Nucl Med. 2020 Apr;61(4):604-607)
% initial values Hill function
p0_hill        = [50 -1];      % see LCN_calc_intact_tracer_hill for details
% boundaries for the models with initial conditions
k0_2T4k_lower_bound        = [0.001 0.001 0.001 0.001];
k0_2T4k_upper_bound        = [1     1     1     1]; 
k0_2T4k_Vb_lower_bound     = [0.001 0.001 0.001 0.001 0.001];
k0_2T4k_Vb_upper_bound     = [1     1     1     1     1]; 
k0_2T4k_vasc1k_lower_bound = [0.001 0.001 0.001 0.001 0.001 0.001];
k0_2T4k_vasc1k_upper_bound = [1     1     1     1     1     1]; 
k0_2T4k_vasc2k_lower_bound = [0.001 0.001 0.001 0.001 0.001 0.001 0.001];
k0_2T4k_vasc2k_upper_bound = [1     1     1     1     1     1     1]; 
nr_params_Logan_input = 2;
nr_params_2T4k    = 4;
nr_params_2T4k_Vb = 5;
nr_2T4k_vasc1k    = 6;
nr_2T4k_vasc2k    = 7;

options = optimoptions('fmincon','Display','off');
STEP = 0.01; % step size for the calculation of integrals in min
% CALC_OPTION: parameter determing the way we calculate the output of a
% model.
%   1 - the model is calculated as the integral of the output concentration
%       devided by the frameduration.
%   2 - the model is calculated at the midscantime.
CALC_OPTION = 2;
window_size = 0; % for temporal median filtering
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
nr_models   = size(model_list,1);


%% READ ATLAS PARCELS
%--------------------------------------------------------------------------

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
      try
         parcel_names = atlas.labels;
      catch
         parcel_names = roi_atlas.labels;
      end
%       parcel_names = num2str(parcel_values);
   end
   nr_VOIS = length(parcel_values);
   
   % initialize output
%    varNames = cell(1,nr_VOIS+1);
%    for voi = 1:nr_VOIS
%        index_value = find(parcel_values_atlas == parcel_values(voi));
%        varNames(1,voi+1) = {[num2str(parcel_values(voi)) '-' char(parcel_names(index_value,:))]};
%    end
   
else
    
   fprintf('You need to define an atlas \n');
   return;
   
end


%% INITIALIZE OUTPUT
%--------------------------------------------------------------------------

% PARAMETERS FOR EACH MODEL

sz = [nr_subjects nr_VOIS+1];
varNames  = [{'subject'} (parcel_names)]; % lixin may26 2025
varTypes(1,1)           = {'string'};
varTypes(1,2:nr_VOIS+1) = {'double'};
results_number_voxels_VOI = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames);
for m = 1:nr_models
    modelname = model_list{m,1};
    if strcmp(modelname,'Logan_input')
       results_region_based_Logan_input_DV    = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames);
       results_region_based_Logan_input_error = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames);
    end
    if strcmp(modelname,'2T4k')
       results_region_based_2T4k_K1    = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames);
       results_region_based_2T4k_k2    = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames);
       results_region_based_2T4k_k3    = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames);
       results_region_based_2T4k_k4    = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames);
       results_region_based_2T4k_DV    = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames);
       results_region_based_2T4k_error = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames);
    end
    if strcmp(modelname,'2T4k_Vb')
       results_region_based_2T4k_Vb_K1    = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames);
       results_region_based_2T4k_Vb_k2    = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames);
       results_region_based_2T4k_Vb_k3    = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames);
       results_region_based_2T4k_Vb_k4    = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames);
       results_region_based_2T4k_Vb_Vb    = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames); 
       results_region_based_2T4k_Vb_DV    = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames);
       results_region_based_2T4k_Vb_error = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames);
    end
    if strcmp(modelname,'2T4k_vasc1k')
       results_region_based_2T4k_vasc1k_K1    = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames);
       results_region_based_2T4k_vasc1k_k2    = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames);
       results_region_based_2T4k_vasc1k_k3    = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames);
       results_region_based_2T4k_vasc1k_k4    = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames);
       results_region_based_2T4k_vasc1k_Vb    = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames); 
       results_region_based_2T4k_vasc1k_K1v   = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames); 
       results_region_based_2T4k_vasc1k_DV    = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames);
       results_region_based_2T4k_vasc1k_error = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames);
    end
    if strcmp(modelname,'2T4k_vasc2k')
       results_region_based_2T4k_vasc2k_K1    = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames);
       results_region_based_2T4k_vasc2k_k2    = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames);
       results_region_based_2T4k_vasc2k_k3    = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames);
       results_region_based_2T4k_vasc2k_k4    = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames);
       results_region_based_2T4k_vasc2k_Vb    = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames); 
       results_region_based_2T4k_vasc2k_K1v   = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames); 
       results_region_based_2T4k_vasc2k_k2v   = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames); 
       results_region_based_2T4k_vasc2k_DV    = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames);
       results_region_based_2T4k_vasc2k_error = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames);
    end
end


% FIT STATISTICS FOR ALL MODELS

sz2 = [nr_subjects*4 nr_VOIS+2];
varNames2  = [{'subject'} {'method'} (parcel_names)]; % lukasvo sept 22 2025
varTypes2 = cell(1, nr_VOIS + 2);
varTypes2(1:2) = {'string', 'string'}; 
varTypes2(3:end) = repmat({'double'}, 1, nr_VOIS);
results_region_based_AIC = table('Size', sz2, 'VariableTypes', varTypes2, 'VariableNames', varNames2);
results_region_based_SC  = table('Size', sz2, 'VariableTypes', varTypes2, 'VariableNames', varNames2);


%% DETERMINING THE GRID OF STARTING VALUES WITHIN THE SEARCH SPACE
%--------------------------------------------------------------------------

for m = 1:nr_models
    modelname = model_list{m,1};
    if strcmp(modelname,'2T4k')
       k0_2T4k_all_values = zeros(nr_starting_values_per_dimension^4,4);
       counter = 0;
       for d1 = 1:nr_starting_values_per_dimension
           val1 = k0_2T4k_lower_bound(1) + d1*(k0_2T4k_upper_bound(1) - k0_2T4k_lower_bound(1))/(nr_starting_values_per_dimension+1);
           for d2 = 1:nr_starting_values_per_dimension
               val2 = k0_2T4k_lower_bound(2) + d2*(k0_2T4k_upper_bound(2) - k0_2T4k_lower_bound(2))/(nr_starting_values_per_dimension+1);
               for d3 = 1:nr_starting_values_per_dimension
                   val3 = k0_2T4k_lower_bound(3) + d3*(k0_2T4k_upper_bound(3) - k0_2T4k_lower_bound(3))/(nr_starting_values_per_dimension+1);
                   for d4 = 1:nr_starting_values_per_dimension
                       val4 = k0_2T4k_lower_bound(4) + d4*(k0_2T4k_upper_bound(4) - k0_2T4k_lower_bound(4))/(nr_starting_values_per_dimension+1);
                       counter = counter + 1;
                       k0_2T4k_all_values(counter,:) = [val1 val2 val3 val4]; 
                   end
               end
           end
       end
       nr_k0_2T4k = size(k0_2T4k_all_values,1);
    end
    if strcmp(modelname,'2T4k_Vb')
       k0_2T4k_Vb_all_values = zeros(nr_starting_values_per_dimension^4,5);
       counter = 0;
       val5 = Vb0;
       for d1 = 1:nr_starting_values_per_dimension
           val1 = k0_2T4k_Vb_lower_bound(1) + d1*(k0_2T4k_Vb_upper_bound(1) - k0_2T4k_Vb_lower_bound(1))/(nr_starting_values_per_dimension+1);
           for d2 = 1:nr_starting_values_per_dimension
               val2 = k0_2T4k_Vb_lower_bound(2) + d2*(k0_2T4k_Vb_upper_bound(2) - k0_2T4k_Vb_lower_bound(2))/(nr_starting_values_per_dimension+1);
               for d3 = 1:nr_starting_values_per_dimension
                   val3 = k0_2T4k_Vb_lower_bound(3) + d3*(k0_2T4k_Vb_upper_bound(3) - k0_2T4k_Vb_lower_bound(3))/(nr_starting_values_per_dimension+1);
                   for d4 = 1:nr_starting_values_per_dimension
                       val4 = k0_2T4k_Vb_lower_bound(4) + d4*(k0_2T4k_Vb_upper_bound(4) - k0_2T4k_Vb_lower_bound(4))/(nr_starting_values_per_dimension+1);
                       counter = counter + 1;
                       k0_2T4k_Vb_all_values(counter,:) = [val1 val2 val3 val4 val5]; 
                   end
               end
           end
       end
       nr_k0_2T4k_Vb = size(k0_2T4k_Vb_all_values,1);
    end
    if strcmp(modelname,'2T4k_vasc1k')
       k0_2T4k_vasc1k_all_values = zeros(nr_starting_values_per_dimension^5,6);
       counter = 0;
       val5 = Vb0;
       for d1 = 1:nr_starting_values_per_dimension
           val1 = k0_2T4k_vasc1k_lower_bound(1) + d1*(k0_2T4k_vasc1k_upper_bound(1) - k0_2T4k_vasc1k_lower_bound(1))/(nr_starting_values_per_dimension+1);
           for d2 = 1:nr_starting_values_per_dimension
               val2 = k0_2T4k_vasc1k_lower_bound(2) + d2*(k0_2T4k_vasc1k_upper_bound(2) - k0_2T4k_vasc1k_lower_bound(2))/(nr_starting_values_per_dimension+1);
               for d3 = 1:nr_starting_values_per_dimension
                   val3 = k0_2T4k_vasc1k_lower_bound(3) + d3*(k0_2T4k_vasc1k_upper_bound(3) - k0_2T4k_vasc1k_lower_bound(3))/(nr_starting_values_per_dimension+1);
                   for d4 = 1:nr_starting_values_per_dimension
                       val4 = k0_2T4k_vasc1k_lower_bound(4) + d4*(k0_2T4k_vasc1k_upper_bound(4) - k0_2T4k_vasc1k_lower_bound(4))/(nr_starting_values_per_dimension+1);
                       for d6 = 1:nr_starting_values_per_dimension
                           val6 = k0_2T4k_vasc1k_lower_bound(6) + d4*(k0_2T4k_vasc1k_upper_bound(6) - k0_2T4k_vasc1k_lower_bound(6))/(nr_starting_values_per_dimension+1);
                           counter = counter + 1;
                           k0_2T4k_vasc1k_all_values(counter,:) = [val1 val2 val3 val4 val5 val6]; 
                       end
                   end
               end
           end
       end
       nr_k0_2T4k_vasc1k = size(k0_2T4k_vasc1k_all_values,1);
    end    
    if strcmp(modelname,'2T4k_vasc2k')
       k0_2T4k_vasc2k_all_values = zeros(nr_starting_values_per_dimension^6,7);
       counter = 0;
       val5 = Vb0;
       for d1 = 1:nr_starting_values_per_dimension
           val1 = k0_2T4k_vasc2k_lower_bound(1) + d1*(k0_2T4k_vasc2k_upper_bound(1) - k0_2T4k_vasc2k_lower_bound(1))/(nr_starting_values_per_dimension+1);
           for d2 = 1:nr_starting_values_per_dimension
               val2 = k0_2T4k_vasc2k_lower_bound(2) + d2*(k0_2T4k_vasc2k_upper_bound(2) - k0_2T4k_vasc2k_lower_bound(2))/(nr_starting_values_per_dimension+1);
               for d3 = 1:nr_starting_values_per_dimension
                   val3 = k0_2T4k_vasc2k_lower_bound(3) + d3*(k0_2T4k_vasc2k_upper_bound(3) - k0_2T4k_vasc2k_lower_bound(3))/(nr_starting_values_per_dimension+1);
                   for d4 = 1:nr_starting_values_per_dimension
                       val4 = k0_2T4k_vasc2k_lower_bound(4) + d4*(k0_2T4k_vasc2k_upper_bound(4) - k0_2T4k_vasc2k_lower_bound(4))/(nr_starting_values_per_dimension+1);
                       for d6 = 1:nr_starting_values_per_dimension
                           val6 = k0_2T4k_vasc2k_lower_bound(6) + d4*(k0_2T4k_vasc2k_upper_bound(6) - k0_2T4k_vasc2k_lower_bound(6))/(nr_starting_values_per_dimension+1);
                           for d7 = 1:nr_starting_values_per_dimension
                               val7 = k0_2T4k_vasc2k_lower_bound(7) + d4*(k0_2T4k_vasc2k_upper_bound(7) - k0_2T4k_vasc2k_lower_bound(7))/(nr_starting_values_per_dimension+1);
                               counter = counter + 1;
                               k0_2T4k_vasc2k_all_values(counter,:) = [val1 val2 val3 val4 val5 val6 val7]; 
                           end
                       end
                   end
               end
           end
       end
       nr_k0_2T4k_vasc2k = size(k0_2T4k_vasc2k_all_values,1);
    end    
end


%% LOOP OVER SUBJECTS
%--------------------------------------------------------------------------

parpool(nr_parpool);

for subj = 1:nr_subjects
    
    % PREP WORK
    
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
    fprintf(fid,'LCN12_PET_TSPO_DPA714_compare_models.m \n');
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
    fprintf(fid,'p0_hill = [%4.2f %4.2f] \n',p0_hill(1),p0_hill(2));

    for m = 1:nr_models
        modelname = model_list{m,1};
        if strcmp(modelname,'Logan_input')
           fprintf(fid,'logan_start_time (min) = %i \n',logan_start_time);
        end
        if strcmp(modelname,'2T4k')
           fprintf(fid,'k0_2T4k_lower_bound = [%4.3f %4.3f %4.3f %4.3f %4.3f] \n',k0_2T4k_lower_bound(1),k0_2T4k_lower_bound(2),k0_2T4k_lower_bound(3),k0_2T4k_lower_bound(4));
           fprintf(fid,'k0_2T4k_upper_bound = [%4.3f %4.3f %4.3f %4.3f %4.3f] \n',k0_2T4k_upper_bound(1),k0_2T4k_upper_bound(2),k0_2T4k_upper_bound(3),k0_2T4k_upper_bound(4));
        end
        if strcmp(modelname,'2T4k_Vb')
           fprintf(fid,'k0_2T4k_Vb_lower_bound = [%4.3f %4.3f %4.3f %4.3f %4.3f] \n',k0_2T4k_Vb_lower_bound(1),k0_2T4k_Vb_lower_bound(2),k0_2T4k_Vb_lower_bound(3),k0_2T4k_Vb_lower_bound(4),k0_2T4k_Vb_lower_bound(5));
           fprintf(fid,'k0_2T4k_Vb_upper_bound = [%4.3f %4.3f %4.3f %4.3f %4.3f] \n',k0_2T4k_Vb_upper_bound(1),k0_2T4k_Vb_upper_bound(2),k0_2T4k_Vb_upper_bound(3),k0_2T4k_Vb_upper_bound(4),k0_2T4k_Vb_upper_bound(5));
        end
        if strcmp(modelname,'2T4k_vasc1k')
           fprintf(fid,'k0_2T4k_vasc1k_lower_bound = [%4.3f %4.3f %4.3f %4.3f %4.3f %4.3f] \n',k0_2T4k_vasc1k_lower_bound(1),k0_2T4k_vasc1k_lower_bound(2),k0_2T4k_vasc1k_lower_bound(3),k0_2T4k_vasc1k_lower_bound(4),k0_2T4k_vasc1k_lower_bound(5),k0_2T4k_vasc1k_lower_bound(6));
           fprintf(fid,'k0_2T4k_vasc1k_upper_bound = [%4.3f %4.3f %4.3f %4.3f %4.3f %4.3f] \n',k0_2T4k_vasc1k_upper_bound(1),k0_2T4k_vasc1k_upper_bound(2),k0_2T4k_vasc1k_upper_bound(3),k0_2T4k_vasc1k_upper_bound(4),k0_2T4k_vasc1k_upper_bound(5),k0_2T4k_vasc1k_upper_bound(6));
        end
        if strcmp(modelname,'2T4k_vasc2k')
           fprintf(fid,'k0_2T4k_vasc2k_lower_bound = [%4.3f %4.3f %4.3f %4.3f %4.3f %4.3f %4.3f] \n',k0_2T4k_vasc2k_lower_bound(1),k0_2T4k_vasc2k_lower_bound(2),k0_2T4k_vasc2k_lower_bound(3),k0_2T4k_vasc2k_lower_bound(4),k0_2T4k_vasc2k_lower_bound(5),k0_2T4k_vasc2k_lower_bound(6),k0_2T4k_vasc2k_lower_bound(7));
           fprintf(fid,'k0_2T4k_vasc2k_upper_bound = [%4.3f %4.3f %4.3f %4.3f %4.3f %4.3f %4.3f] \n',k0_2T4k_vasc2k_upper_bound(1),k0_2T4k_vasc2k_upper_bound(2),k0_2T4k_vasc2k_upper_bound(3),k0_2T4k_vasc2k_upper_bound(4),k0_2T4k_vasc2k_upper_bound(5),k0_2T4k_vasc2k_upper_bound(6),k0_2T4k_vasc2k_upper_bound(7));
        end
    end

    
    % CHECK IF WE CAN FIND ALL INPUT DATA IN MNI SPACE
    
    [filename_GM,go_GM]    = LCN_check_filename(dir_deriv_anat2,['wc1' subjectname '*.nii']);
    [filename_WM,go_WM]    = LCN_check_filename(dir_deriv_anat2,['wc2' subjectname '*.nii']);
    [filename_CSF,go_CSF]  = LCN_check_filename(dir_deriv_anat2,['wc3' subjectname '*.nii']);
    [filename_wPET,go_PET] = LCN_check_filename(dir_deriv,['w' subjectname '*_' infostring_tracer '*_' infostring_PET '.nii']); 
    [frame_definition_file,go_frames] = LCN_check_filename(dir_pet,[subjectname '_*' infostring_tracer '_' infostring_PET '_' infostring_frames '.m']); 
    [input_file,go_input]  = LCN_check_filename(dir_pet,[subjectname '_*' infostring_tracer '_' infostring_input '.m']); 
    [metab_file,go_metab]  = LCN_check_filename(dir_pet,[subjectname '_*' infostring_tracer '_' infostring_metab '.m']); 

    
    % READ PLASMA AND BLOOD DATA, PERFORM CALCULATIONS, AND ANALYSE DATA
    
    if go_frames == 1 % frame definition file found   
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
    
    
    % READ BRAIN DATA
    
    if go_PET == 1 % dynamic data are found   
       % reading data in MNI space
       fprintf('reading data of subject %s',subjectname);
       fprintf(fid,'\n');
       fprintf(fid,'Reading data in MNI space started: %s \n',datetime('now'));
       
       if pet_space_reference
        Vref_tmp = spm_vol(filename_wPET); 
        Vref = Vref_tmp(1);
        % read dynamic data
        dydata = zeros(Vref.dim(1),Vref.dim(2),Vref.dim(3),nr_frames);
        for frame = 1:nr_frames
           clear tmp
           tmp = LCN12_read_image([filename_wPET ',' num2str(frame)],Vref);
           dydata(:,:,:,frame) = tmp/1000; %in kBq/ml
        end
       
       else
         atlas = atlas_orig;
         Vref = Vatlas_orig;
         dydata = zeros(Vref.dim(1),Vref.dim(2),Vref.dim(3),nr_frames);
         for frame = 1:nr_frames
            clear tmp
            tmp = LCN12_read_image([filename_wPET ',' num2str(frame)],Vref);
            dydata(:,:,:,frame) = tmp/1000; %in kBq/ml
         end
         
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
    
    
    % INITIATE VARIABLES FOR MODEL (COMPARISON) RESULTS
    
    for m = 1:nr_models
        modelname = model_list{m,1};
        if strcmp(modelname,'Logan_input')
           results_region_based_Logan_input_DV(subj,1)    = {subjectname};
           results_region_based_Logan_input_error(subj,1) = {subjectname};
           if save_parcelimgs
               atlas_Logan_input_DV    = zeros(size(atlas_orig));
               atlas_Logan_input_error = zeros(size(atlas_orig));
           end
        end
        if strcmp(modelname,'2T4k')
           results_region_based_2T4k_K1(subj,1)    = {subjectname};
           results_region_based_2T4k_k2(subj,1)    = {subjectname};
           results_region_based_2T4k_k3(subj,1)    = {subjectname};
           results_region_based_2T4k_k4(subj,1)    = {subjectname};
           results_region_based_2T4k_DV(subj,1)    = {subjectname};
           results_region_based_2T4k_error(subj,1) = {subjectname};
           if save_parcelimgs
               atlas_2T4k_K1    = zeros(size(atlas_orig));
               atlas_2T4k_k2    = zeros(size(atlas_orig));
               atlas_2T4k_k3    = zeros(size(atlas_orig));
               atlas_2T4k_k4    = zeros(size(atlas_orig));
               atlas_2T4k_DV    = zeros(size(atlas_orig));
               atlas_2T4k_error = zeros(size(atlas_orig));
           end
        end
        if strcmp(modelname,'2T4k_Vb')
           results_region_based_2T4k_Vb_K1(subj,1)    = {subjectname};
           results_region_based_2T4k_Vb_k2(subj,1)    = {subjectname};
           results_region_based_2T4k_Vb_k3(subj,1)    = {subjectname};
           results_region_based_2T4k_Vb_k4(subj,1)    = {subjectname};
           results_region_based_2T4k_Vb_Vb(subj,1)    = {subjectname};
           results_region_based_2T4k_Vb_DV(subj,1)    = {subjectname};
           results_region_based_2T4k_Vb_error(subj,1) = {subjectname};
           if save_parcelimgs
               atlas_2T4k_Vb_K1    = zeros(size(atlas_orig));
               atlas_2T4k_Vb_k2    = zeros(size(atlas_orig));
               atlas_2T4k_Vb_k3    = zeros(size(atlas_orig));
               atlas_2T4k_Vb_k4    = zeros(size(atlas_orig));
               atlas_2T4k_Vb_Vb    = zeros(size(atlas_orig));
               atlas_2T4k_Vb_DV    = zeros(size(atlas_orig));
               atlas_2T4k_Vb_error = zeros(size(atlas_orig));
           end
        end
        if strcmp(modelname,'2T4k_vasc1k')
           results_region_based_2T4k_vasc1k_K1(subj,1)    = {subjectname};
           results_region_based_2T4k_vasc1k_k2(subj,1)    = {subjectname};
           results_region_based_2T4k_vasc1k_k3(subj,1)    = {subjectname};
           results_region_based_2T4k_vasc1k_k4(subj,1)    = {subjectname};
           results_region_based_2T4k_vasc1k_Vb(subj,1)    = {subjectname};
           results_region_based_2T4k_vasc1k_K1v(subj,1)   = {subjectname};
           results_region_based_2T4k_vasc1k_DV(subj,1)    = {subjectname};
           results_region_based_2T4k_vasc1k_error(subj,1) = {subjectname};
           if save_parcelimgs
               atlas_2T4k_vasc1k_K1    = zeros(size(atlas_orig));
               atlas_2T4k_vasc1k_k2    = zeros(size(atlas_orig));
               atlas_2T4k_vasc1k_k3    = zeros(size(atlas_orig));
               atlas_2T4k_vasc1k_k4    = zeros(size(atlas_orig));
               atlas_2T4k_vasc1k_Vb    = zeros(size(atlas_orig));
               atlas_2T4k_vasc1k_K1v   = zeros(size(atlas_orig));
               atlas_2T4k_vasc1k_DV    = zeros(size(atlas_orig));
               atlas_2T4k_vasc1k_error = zeros(size(atlas_orig));
           end
        end
        if strcmp(modelname,'2T4k_vasc2k')
           results_region_based_2T4k_vasc2k_K1(subj,1)    = {subjectname};
           results_region_based_2T4k_vasc2k_k2(subj,1)    = {subjectname};
           results_region_based_2T4k_vasc2k_k3(subj,1)    = {subjectname};
           results_region_based_2T4k_vasc2k_k4(subj,1)    = {subjectname};
           results_region_based_2T4k_vasc2k_Vb(subj,1)    = {subjectname};
           results_region_based_2T4k_vasc2k_K1v(subj,1)   = {subjectname};
           results_region_based_2T4k_vasc2k_k2v(subj,1)   = {subjectname};
           results_region_based_2T4k_vasc2k_DV(subj,1)    = {subjectname};
           results_region_based_2T4k_vasc2k_error(subj,1) = {subjectname};
           if save_parcelimgs
               atlas_2T4k_vasc2k_K1    = zeros(size(atlas_orig));
               atlas_2T4k_vasc2k_k2    = zeros(size(atlas_orig));
               atlas_2T4k_vasc2k_k3    = zeros(size(atlas_orig));
               atlas_2T4k_vasc2k_k4    = zeros(size(atlas_orig));
               atlas_2T4k_vasc2k_Vb    = zeros(size(atlas_orig));
               atlas_2T4k_vasc2k_K1v   = zeros(size(atlas_orig));
               atlas_2T4k_vasc2k_k2v   = zeros(size(atlas_orig));
               atlas_2T4k_vasc2k_DV    = zeros(size(atlas_orig));
               atlas_2T4k_vasc2k_error = zeros(size(atlas_orig));
           end
        end
    end
    
    results_number_voxels_VOI(subj,1) = {subjectname};
    results_region_based_AIC(4*(subj-1)+1,1)  = {subjectname};
    results_region_based_SC(4*(subj-1)+1,1)   = {subjectname};

    
    % READ ATLAS DATA
    
    if pet_space_reference == 1
    % read the atlas in the matrix of the subject
    atlas = LCN12_read_image(atlas_filename,Vref);
    % atlas contains non-integer values in voxels which are a mixture. We focus
    % only on voxels which are of one specific type
    % therefore, set the value of other voxels to 0
    atlas(mod(atlas,1) > 0) = 0;
    end
    % adapt weight of the frames 
    WEIGHTS_PET = ones(length(MIDSCANTIMES),1);
    WEIGHTS_PET(excluded_frames) = 0;
    WEIGHTS_PET    = WEIGHTS_PET./sum(WEIGHTS_PET);
    fprintf('atlas based analysis started at %s \n',datetime('now'));
    fprintf(fid,'atlas based analysis started at %s \n',datetime('now'));
    
    
    % LOOP OVER VOIS
    
    for voi = 1:nr_VOIS
        clear index_value VOIimg VOIname C_MEASURED
        index_value = find(parcel_values_atlas == parcel_values(voi));
        VOIname  = parcel_names{1,index_value};
        fprintf(fid,'VOI %i: %s\n',index_value,VOIname);
        VOIimg = (atlas == parcel_values(voi));
        % intersect VOI with subject specific GM or WM
        if intersect_GM == 1
           VOIimg = (VOIimg > 0.5).*(GMimg > GM_CUTOFF);
        elseif intersect_WM  == 1
           VOIimg = (VOIimg > 0.5).*(WMimg > WM_CUTOFF);
        end
        final_VOI_mask = (VOIimg > 0.5).*brain_mask > 0;
        results_number_voxels_VOI(subj,1+voi) = {nnz(final_VOI_mask(:) > 0)};
        if nnz(final_VOI_mask(:) > 0) < min_number_voxels 
           continue;
        end
        C_MEASURED = zeros(nr_frames,1);
        for frame = 1:nr_frames
            clear tmp values
            tmp = dydata(:,:,:,frame);
            values = tmp(final_VOI_mask > 0);
            C_MEASURED(frame) = mean(values(~isnan(values)));
        end
        if window_size > 0
           % temporal filtering using a median filter
           C_MEASURED = medfilt1(C_MEASURED, window_size);
        end
        if sum(abs(C_MEASURED)) == 0
           continue;
        end
        if sum(isnan(C_MEASURED)) > 0
           continue;
        end
        
        
        % CALCULATE ALL THE MODELS
        
        if figures_on == 1 || save_figures == 1 
           close(3);
           hfig3 = figure(3);
           set(hfig3,'Name',[subjectname ' - ' VOIname]);
           close(4);
           hfig4 = figure(4);
           set(hfig4,'Name',[subjectname ' - ' VOIname]);      
        end
        
        for m = 1:nr_models
            modelname = model_list{m,1};
            if strcmp(modelname,'Logan_input')
               clear VD error
               if figures_on == 1 || save_figures == 1
                  figure(3);
                  [DV,error] = LCN_LOGAN(MIDSCANTIMES(frames_ok),C_MEASURED(frames_ok),time_input,corrected_input,logan_start_time,3);
                  title('Logan input');
                  if save_figures == 1
                     cd(dir_figures);
                     saveas(hfig3,['fig_Logan_input_VOI_' VOIname '_' suffix '.fig']);
                  else 
                     pause
                  end
               else
                  [DV,error] = LCN_LOGAN(MIDSCANTIMES(frames_ok),C_MEASURED(frames_ok),time_input,corrected_input,logan_start_time,0);
               end
               results_region_based_Logan_input_DV(subj,1+voi)    = {DV};
               results_region_based_Logan_input_error(subj,1+voi) = {error};
               if save_parcelimgs
                   atlas_Logan_input_DV(atlas_orig == parcel_values(voi))    = DV;
                   atlas_Logan_input_error(atlas_orig == parcel_values(voi)) = error;
                   Logan_outputname1 = fullfile(dir_firstlevel,['w' subjectname '_' suffix '_atlas_Logan_input_DV.nii']);
                   Logan_outputname2 = fullfile(dir_firstlevel,['w' subjectname '_' suffix '_atlas_Logan_input_error.nii']);
                   LCN12_write_image(atlas_Logan_input_DV,Logan_outputname1, '_' + suffix + '_Logan input DV',64,Vatlas_orig);        
                   LCN12_write_image(atlas_Logan_input_error,Logan_outputname2, '_' + suffix + '_Logan input error',64,Vatlas_orig);
               end
            end
        
            if strcmp(modelname,'2T4k')
               clear tmpA tmpB index_min kfit_2T4k cost_function axvalues k0
               % start new fit from different starting positions
               tmpA = zeros(nr_k0_2T4k,4);
               tmpB = zeros(nr_k0_2T4k,1);
               cost_function = @(params) norm(WEIGHTS_PET.*(LCN_calc2_model_2T4k(params,MIDSCANTIMES,FRAMEDURATION,CA_FINETIMES,FINETIMES,CALC_OPTION,STEP) - C_MEASURED));
               parfor j = 1:nr_k0_2T4k
                      k0 = k0_2T4k_all_values(j,:);
                      [tmpA(j,:),tmpB(j)] = fmincon(cost_function,k0,[],[],[],[],k0_2T4k_lower_bound,k0_2T4k_upper_bound,[],options);
               end
               index_min = find(tmpB == min(tmpB));
               kfit_2T4k = tmpA(index_min(1),:);
               results_region_based_2T4k_K1(subj,1+voi)    = {kfit_2T4k(1)};
               results_region_based_2T4k_k2(subj,1+voi)    = {kfit_2T4k(2)};
               results_region_based_2T4k_k3(subj,1+voi)    = {kfit_2T4k(3)};
               results_region_based_2T4k_k4(subj,1+voi)    = {kfit_2T4k(4)};
               results_region_based_2T4k_DV(subj,1+voi)    = {kfit_2T4k(1)*(1+kfit_2T4k(3)/kfit_2T4k(4))/kfit_2T4k(2)};
               results_region_based_2T4k_error(subj,1+voi) = {tmpB(index_min(1))./mean(C_MEASURED)}; % relative error
               if save_parcelimgs
                   atlas_2T4k_K1(atlas_orig == parcel_values(voi))    = kfit_2T4k(1);
                   atlas_2T4k_k2(atlas_orig == parcel_values(voi))    = kfit_2T4k(2);
                   atlas_2T4k_k3(atlas_orig == parcel_values(voi))    = kfit_2T4k(3);
                   atlas_2T4k_k4(atlas_orig == parcel_values(voi))    = kfit_2T4k(4);
                   atlas_2T4k_DV(atlas_orig == parcel_values(voi))    = kfit_2T4k(1)*(1+kfit_2T4k(3)/kfit_2T4k(4))/kfit_2T4k(2);
                   atlas_2T4k_error(atlas_orig == parcel_values(voi)) = tmpB(index_min(1))./mean(C_MEASURED);
                   twoT4k_outputname1 = fullfile(dir_firstlevel,['w' subjectname '_' suffix '_atlas_2T4k_K1.nii']);
                   twoT4k_outputname2 = fullfile(dir_firstlevel,['w' subjectname '_' suffix '_atlas_2T4k_k2.nii']);
                   twoT4k_outputname3 = fullfile(dir_firstlevel,['w' subjectname '_' suffix '_atlas_2T4k_k3.nii']);
                   twoT4k_outputname4 = fullfile(dir_firstlevel,['w' subjectname '_' suffix '_atlas_2T4k_k4.nii']);
                   twoT4k_outputname6 = fullfile(dir_firstlevel,['w' subjectname '_' suffix '_atlas_2T4k_DV.nii']);
                   twoT4k_outputname7 = fullfile(dir_firstlevel,['w' subjectname '_' suffix '_atlas_2T4k_error.nii']);
                   LCN12_write_image(atlas_2T4k_K1,twoT4k_outputname1,'_' + suffix + '_2T4k_K1',64,Vatlas_orig);        
                   LCN12_write_image(atlas_2T4k_k2,twoT4k_outputname2,'_' + suffix + '_2T4k_k2',64,Vatlas_orig);        
                   LCN12_write_image(atlas_2T4k_k3,twoT4k_outputname3,'_' + suffix + '_2T4k_k3',64,Vatlas_orig);        
                   LCN12_write_image(atlas_2T4k_k4,twoT4k_outputname4,'_' + suffix + '_2T4k_k4',64,Vatlas_orig);               
                   LCN12_write_image(atlas_2T4k_DV,twoT4k_outputname6,'_' + suffix + '_2T4k_DV',64,Vatlas_orig);        
                   LCN12_write_image(atlas_2T4k_error,twoT4k_outputname7,'_' + suffix + '_2T4k_error',64,Vatlas_orig);  
               end
               [AIC,SC] = LCN_calc_model_selection(nnz(WEIGHTS_PET),tmpB(index_min(1)).^2,nr_params_2T4k);
               results_region_based_AIC((subj-1)*4+1,2) = {modelname};
               results_region_based_SC((subj-1)*4+1,2)  = {modelname};
               results_region_based_AIC((subj-1)*4+1,2+voi) = {AIC};
               results_region_based_SC((subj-1)*4+1,2+voi)  = {SC};
               if figures_on == 1 || save_figures == 1
                  figure(4)
                  subplot(221)
                  plot(MIDSCANTIMES,C_MEASURED,'ob:');
                  hold on
                  plot(MIDSCANTIMES(excluded_frames),C_MEASURED(excluded_frames),'*r'); % indicating unreliable data due to head movements
                  xlabel('time in min');
                  ylabel('kBq/ml');
                  title('2T4k');
                  plot(MIDSCANTIMES,LCN_calc2_model_2T4k(kfit_2T4k,MIDSCANTIMES,FRAMEDURATION,CA_FINETIMES,FINETIMES,CALC_OPTION,STEP),'k')
                  axvalues = axis;
                  text(axvalues(2)*0.5,axvalues(3)+0.1*(axvalues(4)-axvalues(3)),['VD = ' num2str(kfit_2T4k(1)*(1+kfit_2T4k(3)/kfit_2T4k(4))/kfit_2T4k(2))]);
               end
            end
            
            if strcmp(modelname,'2T4k_Vb')
               clear tmpA tmpB index_min kfit_2T4k_Vb cost_function axvalues k0
               % start new fit from different starting positions
               tmpA = zeros(nr_k0_2T4k_Vb,5);
               tmpB = zeros(nr_k0_2T4k_Vb,1);
               cost_function = @(params) norm(WEIGHTS_PET.*(LCN_calc2_model_2T4k_Vb(params,MIDSCANTIMES,FRAMEDURATION,CA_FINETIMES,CA_WB_FINETIMES,FINETIMES,CALC_OPTION,STEP) - C_MEASURED));
               parfor j = 1:nr_k0_2T4k_Vb
                      k0 = k0_2T4k_Vb_all_values(j,:);
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
               if save_parcelimgs
                   atlas_2T4k_Vb_K1(atlas_orig == parcel_values(voi))    = kfit_2T4k_Vb(1);
                   atlas_2T4k_Vb_k2(atlas_orig == parcel_values(voi))    = kfit_2T4k_Vb(2);
                   atlas_2T4k_Vb_k3(atlas_orig == parcel_values(voi))    = kfit_2T4k_Vb(3);
                   atlas_2T4k_Vb_k4(atlas_orig == parcel_values(voi))    = kfit_2T4k_Vb(4);
                   atlas_2T4k_Vb_Vb(atlas_orig == parcel_values(voi))    = kfit_2T4k_Vb(5);
                   atlas_2T4k_Vb_DV(atlas_orig == parcel_values(voi))    = kfit_2T4k_Vb(1)*(1+kfit_2T4k_Vb(3)/kfit_2T4k_Vb(4))/kfit_2T4k_Vb(2);
                   atlas_2T4k_Vb_error(atlas_orig == parcel_values(voi)) = tmpB(index_min(1))./mean(C_MEASURED);
                   twoT4k_Vb_outputname1 = fullfile(dir_firstlevel,['w' subjectname '_' suffix '_atlas_2T4k_Vb_K1.nii']);
                   twoT4k_Vb_outputname2 = fullfile(dir_firstlevel,['w' subjectname '_' suffix '_atlas_2T4k_Vb_k2.nii']);
                   twoT4k_Vb_outputname3 = fullfile(dir_firstlevel,['w' subjectname '_' suffix '_atlas_2T4k_Vb_k3.nii']);
                   twoT4k_Vb_outputname4 = fullfile(dir_firstlevel,['w' subjectname '_' suffix '_atlas_2T4k_Vb_k4.nii']);
                   twoT4k_Vb_outputname5 = fullfile(dir_firstlevel,['w' subjectname '_' suffix '_atlas_2T4k_Vb_Vb.nii']);
                   twoT4k_Vb_outputname6 = fullfile(dir_firstlevel,['w' subjectname '_' suffix '_atlas_2T4k_Vb_DV.nii']);
                   twoT4k_Vb_outputname7 = fullfile(dir_firstlevel,['w' subjectname '_' suffix '_atlas_2T4k_Vb_error.nii']);
                   LCN12_write_image(atlas_2T4k_Vb_K1,twoT4k_Vb_outputname1,'_' + suffix + '2T4k_Vb_K1',64,Vatlas_orig);        
                   LCN12_write_image(atlas_2T4k_Vb_k2,twoT4k_Vb_outputname2,'_' + suffix + '2T4k_Vb_k2',64,Vatlas_orig);        
                   LCN12_write_image(atlas_2T4k_Vb_k3,twoT4k_Vb_outputname3,'_' + suffix + '2T4k_Vb_k3',64,Vatlas_orig);        
                   LCN12_write_image(atlas_2T4k_Vb_k4,twoT4k_Vb_outputname4,'_' + suffix + '2T4k_Vb_k4',64,Vatlas_orig);        
                   LCN12_write_image(atlas_2T4k_Vb_Vb,twoT4k_Vb_outputname5,'_' + suffix + '2T4k_Vb_Vb',64,Vatlas_orig);        
                   LCN12_write_image(atlas_2T4k_Vb_DV,twoT4k_Vb_outputname6,'_' + suffix + '2T4k_Vb_DV',64,Vatlas_orig);        
                   LCN12_write_image(atlas_2T4k_Vb_error,twoT4k_Vb_outputname7,'_' + suffix + '2T4k_Vb_error',64,Vatlas_orig); 
               end
               [AIC,SC] = LCN_calc_model_selection(nnz(WEIGHTS_PET),tmpB(index_min(1)).^2,nr_params_2T4k_Vb);
               results_region_based_AIC((subj-1)*4+2,2) = {modelname};
               results_region_based_SC((subj-1)*4+2,2)  = {modelname};
               results_region_based_AIC((subj-1)*4+2,2+voi) = {AIC};
               results_region_based_SC((subj-1)*4+2,2+voi)  = {SC};
               if figures_on == 1 || save_figures == 1 
                  figure(4)
                  subplot(222)
                  plot(MIDSCANTIMES,C_MEASURED,'ob:');
                  hold on
                  plot(MIDSCANTIMES(excluded_frames),C_MEASURED(excluded_frames),'*r'); % indicating unreliable data due to head movements
                  xlabel('time in min');
                  ylabel('kBq/ml');
                  title('2T4k\_Vb');
                  plot(MIDSCANTIMES,LCN_calc2_model_2T4k_Vb(kfit_2T4k_Vb,MIDSCANTIMES,FRAMEDURATION,CA_FINETIMES,CA_WB_FINETIMES,FINETIMES,CALC_OPTION,STEP),'k')
                  axvalues = axis;
                  text(axvalues(2)*0.5,axvalues(3)+0.1*(axvalues(4)-axvalues(3)),['VD = ' num2str(kfit_2T4k_Vb(1)*(1+kfit_2T4k_Vb(3)/kfit_2T4k_Vb(4))/kfit_2T4k_Vb(2))]);
               end
            end
            
            if strcmp(modelname,'2T4k_vasc1k')
               clear tmpA tmpB index_min kfit_2T4k_vasc1k cost_function axvalues k0
               % start new fit from different starting positions
               tmpA = zeros(nr_k0_2T4k_vasc1k,6);
               tmpB = zeros(nr_k0_2T4k_vasc1k,1);
               cost_function = @(params) norm(WEIGHTS_PET.*(LCN_calc_model_2T4k_vasc1k(params,MIDSCANTIMES,FRAMEDURATION,CA_FINETIMES,CA_WB_FINETIMES,FINETIMES,CALC_OPTION,STEP) - C_MEASURED));
               parfor j = 1:nr_k0_2T4k_vasc1k
                      k0 = k0_2T4k_vasc1k_all_values(j,:);
                      [tmpA(j,:),tmpB(j)] = fmincon(cost_function,k0,[],[],[],[],k0_2T4k_vasc1k_lower_bound,k0_2T4k_vasc1k_upper_bound,[],options);
               end
               index_min = find(tmpB == min(tmpB));
               kfit_2T4k_vasc1k = tmpA(index_min(1),:);
               results_region_based_2T4k_vasc1k_K1(subj,1+voi)    = {kfit_2T4k_vasc1k(1)};
               results_region_based_2T4k_vasc1k_k2(subj,1+voi)    = {kfit_2T4k_vasc1k(2)};
               results_region_based_2T4k_vasc1k_k3(subj,1+voi)    = {kfit_2T4k_vasc1k(3)};
               results_region_based_2T4k_vasc1k_k4(subj,1+voi)    = {kfit_2T4k_vasc1k(4)};
               results_region_based_2T4k_vasc1k_Vb(subj,1+voi)    = {kfit_2T4k_vasc1k(5)};
               results_region_based_2T4k_vasc1k_K1v(subj,1+voi)   = {kfit_2T4k_vasc1k(6)};
               results_region_based_2T4k_vasc1k_DV(subj,1+voi)    = {kfit_2T4k_vasc1k(1)*(1+kfit_2T4k_vasc1k(3)/kfit_2T4k_vasc1k(4))/kfit_2T4k_vasc1k(2)};
               results_region_based_2T4k_vasc1k_error(subj,1+voi) = {tmpB(index_min(1))./mean(C_MEASURED)}; % relative error
               if save_parcelimgs
                   atlas_2T4k_vasc1k_K1(atlas_orig == parcel_values(voi))    = kfit_2T4k_vasc1k(1);
                   atlas_2T4k_vasc1k_k2(atlas_orig == parcel_values(voi))    = kfit_2T4k_vasc1k(2);
                   atlas_2T4k_vasc1k_k3(atlas_orig == parcel_values(voi))    = kfit_2T4k_vasc1k(3);
                   atlas_2T4k_vasc1k_k4(atlas_orig == parcel_values(voi))    = kfit_2T4k_vasc1k(4);
                   atlas_2T4k_vasc1k_Vb(atlas_orig == parcel_values(voi))    = kfit_2T4k_vasc1k(5);
                   atlas_2T4k_vasc1k_K1v(atlas_orig == parcel_values(voi))   = kfit_2T4k_vasc1k(6);
                   atlas_2T4k_vasc1k_DV(atlas_orig == parcel_values(voi))    = kfit_2T4k_vasc1k(1)*(1+kfit_2T4k_vasc1k(3)/kfit_2T4k_vasc1k(4))/kfit_2T4k_vasc1k(2);
                   atlas_2T4k_vasc1k_error(atlas_orig == parcel_values(voi)) = tmpB(index_min(1))./mean(C_MEASURED);
                   twoT4k_vasc1k_outputname1 = fullfile(dir_firstlevel,['w' subjectname '_' suffix '_atlas_2T4k_vasc1k_K1.nii']);
                   twoT4k_vasc1k_outputname2 = fullfile(dir_firstlevel,['w' subjectname '_' suffix '_atlas_2T4k_vasc1k_k2.nii']);
                   twoT4k_vasc1k_outputname3 = fullfile(dir_firstlevel,['w' subjectname '_' suffix '_atlas_2T4k_vasc1k_k3.nii']);
                   twoT4k_vasc1k_outputname4 = fullfile(dir_firstlevel,['w' subjectname '_' suffix '_atlas_2T4k_vasc1k_k4.nii']);
                   twoT4k_vasc1k_outputname5 = fullfile(dir_firstlevel,['w' subjectname '_' suffix '_atlas_2T4k_vasc1k_Vb.nii']);
                   twoT4k_vasc1k_outputname6 = fullfile(dir_firstlevel,['w' subjectname '_' suffix '_atlas_2T4k_vasc1k_K1v.nii']);
                   twoT4k_vasc1k_outputname7 = fullfile(dir_firstlevel,['w' subjectname '_' suffix '_atlas_2T4k_vasc1k_DV.nii']);
                   twoT4k_vasc1k_outputname8 = fullfile(dir_firstlevel,['w' subjectname '_' suffix '_atlas_2T4k_vasc1k_error.nii']);
                   LCN12_write_image(atlas_2T4k_vasc1k_K1,twoT4k_vasc1k_outputname1,'_' + suffix + '2T4k_vasc1k_K1',64,Vatlas_orig);        
                   LCN12_write_image(atlas_2T4k_vasc1k_k2,twoT4k_vasc1k_outputname2,'_' + suffix + '2T4k_vasc1k_k2',64,Vatlas_orig);        
                   LCN12_write_image(atlas_2T4k_vasc1k_k3,twoT4k_vasc1k_outputname3,'_' + suffix + '2T4k_vasc1k_k3',64,Vatlas_orig);        
                   LCN12_write_image(atlas_2T4k_vasc1k_k4,twoT4k_vasc1k_outputname4,'_' + suffix + '2T4k_vasc1k_k4',64,Vatlas_orig);        
                   LCN12_write_image(atlas_2T4k_vasc1k_Vb,twoT4k_vasc1k_outputname5,'_' + suffix + '2T4k_vasc1k_Vb',64,Vatlas_orig);  
                   LCN12_write_image(atlas_2T4k_vasc1k_K1v,twoT4k_vasc1k_outputname6,'_' + suffix + '2T4k_vasc1k_K1v',64,Vatlas_orig);
                   LCN12_write_image(atlas_2T4k_vasc1k_DV,twoT4k_vasc1k_outputname7,'_' + suffix + '2T4k_vasc1k_DV',64,Vatlas_orig);        
                   LCN12_write_image(atlas_2T4k_vasc1k_error,twoT4k_vasc1k_outputname8,'_' + suffix + '2T4k_vasc1k_error',64,Vatlas_orig);
               end
               [AIC,SC] = LCN_calc_model_selection(nnz(WEIGHTS_PET),tmpB(index_min(1)).^2,nr_params_2T4k_Vb);
               results_region_based_AIC((subj-1)*4+3,2) = {modelname};
               results_region_based_SC((subj-1)*4+3,2)  = {modelname};
               results_region_based_AIC((subj-1)*4+3,2+voi) = {AIC};
               results_region_based_SC((subj-1)*4+3,2+voi)  = {SC};
               if figures_on == 1 || save_figures == 1 
                  figure(4)
                  subplot(223)
                  plot(MIDSCANTIMES,C_MEASURED,'ob:');
                  hold on
                  plot(MIDSCANTIMES(excluded_frames),C_MEASURED(excluded_frames),'*r'); % indicating unreliable data due to head movements
                  xlabel('time in min');
                  ylabel('kBq/ml');
                  title('2T4K\_vasc1k');
                  plot(MIDSCANTIMES,LCN_calc_model_2T4k_vasc1k(kfit_2T4k_vasc1k,MIDSCANTIMES,FRAMEDURATION,CA_FINETIMES,CA_WB_FINETIMES,FINETIMES,CALC_OPTION,STEP),'k')
                  axvalues = axis;
                  text(axvalues(2)*0.5,axvalues(3)+0.1*(axvalues(4)-axvalues(3)),['VD = ' num2str(kfit_2T4k_vasc1k(1)*(1+kfit_2T4k_vasc1k(3)/kfit_2T4k_vasc1k(4))/kfit_2T4k_vasc1k(2))]);
               end            
            end
            
            if strcmp(modelname,'2T4k_vasc2k')
               clear tmpA tmpB index_min kfit_2T4k_vasc2k cost_function axvalues k0
               % start new fit from different starting positions
               tmpA = zeros(nr_k0_2T4k_vasc2k,7);
               tmpB = zeros(nr_k0_2T4k_vasc2k,1);
               cost_function = @(params) norm(WEIGHTS_PET.*(LCN_calc_model_2T4k_vasc2k(params,MIDSCANTIMES,FRAMEDURATION,CA_FINETIMES,CA_WB_FINETIMES,FINETIMES,CALC_OPTION,STEP) - C_MEASURED));
               parfor j = 1:nr_k0_2T4k_vasc2k
                      k0 = k0_2T4k_vasc2k_all_values(j,:);
                      [tmpA(j,:),tmpB(j)] = fmincon(cost_function,k0,[],[],[],[],k0_2T4k_vasc2k_lower_bound,k0_2T4k_vasc2k_upper_bound,[],options);
               end
               index_min = find(tmpB == min(tmpB));
               kfit_2T4k_vasc2k = tmpA(index_min(1),:);
               results_region_based_2T4k_vasc2k_K1(subj,1+voi)    = {kfit_2T4k_vasc2k(1)};
               results_region_based_2T4k_vasc2k_k2(subj,1+voi)    = {kfit_2T4k_vasc2k(2)};
               results_region_based_2T4k_vasc2k_k3(subj,1+voi)    = {kfit_2T4k_vasc2k(3)};
               results_region_based_2T4k_vasc2k_k4(subj,1+voi)    = {kfit_2T4k_vasc2k(4)};
               results_region_based_2T4k_vasc2k_Vb(subj,1+voi)    = {kfit_2T4k_vasc2k(5)};
               results_region_based_2T4k_vasc2k_K1v(subj,1+voi)   = {kfit_2T4k_vasc2k(6)};
               results_region_based_2T4k_vasc2k_k2v(subj,1+voi)   = {kfit_2T4k_vasc2k(7)};
               results_region_based_2T4k_vasc2k_DV(subj,1+voi)    = {kfit_2T4k_vasc2k(1)*(1+kfit_2T4k_vasc2k(3)/kfit_2T4k_vasc2k(4))/kfit_2T4k_vasc2k(2)};
               results_region_based_2T4k_vasc2k_error(subj,1+voi) = {tmpB(index_min(1))./mean(C_MEASURED)}; % relative error
               if save_parcelimgs
                   atlas_2T4k_vasc2k_K1(atlas_orig == parcel_values(voi))    = kfit_2T4k_vasc2k(1);
                   atlas_2T4k_vasc2k_k2(atlas_orig == parcel_values(voi))    = kfit_2T4k_vasc2k(2);
                   atlas_2T4k_vasc2k_k3(atlas_orig == parcel_values(voi))    = kfit_2T4k_vasc2k(3);
                   atlas_2T4k_vasc2k_k4(atlas_orig == parcel_values(voi))    = kfit_2T4k_vasc2k(4);
                   atlas_2T4k_vasc2k_Vb(atlas_orig == parcel_values(voi))    = kfit_2T4k_vasc2k(5);
                   atlas_2T4k_vasc2k_K1v(atlas_orig == parcel_values(voi))   = kfit_2T4k_vasc2k(6);
                   atlas_2T4k_vasc2k_k2v(atlas_orig == parcel_values(voi))   = kfit_2T4k_vasc2k(7);
                   atlas_2T4k_vasc2k_DV(atlas_orig == parcel_values(voi))    = kfit_2T4k_vasc2k(1)*(1+kfit_2T4k_vasc2k(3)/kfit_2T4k_vasc2k(4))/kfit_2T4k_vasc2k(2);
                   atlas_2T4k_vasc2k_error(atlas_orig == parcel_values(voi)) = tmpB(index_min(1))./mean(C_MEASURED);
                   twoT4k_vasc2k_outputname1 = fullfile(dir_firstlevel,['w' subjectname '_' suffix '_atlas_2T4k_vasc2k_K1.nii']);
                   twoT4k_vasc2k_outputname2 = fullfile(dir_firstlevel,['w' subjectname '_' suffix '_atlas_2T4k_vasc2k_k2.nii']);
                   twoT4k_vasc2k_outputname3 = fullfile(dir_firstlevel,['w' subjectname '_' suffix '_atlas_2T4k_vasc2k_k3.nii']);
                   twoT4k_vasc2k_outputname4 = fullfile(dir_firstlevel,['w' subjectname '_' suffix '_atlas_2T4k_vasc2k_k4.nii']);
                   twoT4k_vasc2k_outputname5 = fullfile(dir_firstlevel,['w' subjectname '_' suffix '_atlas_2T4k_vasc2k_Vb.nii']);
                   twoT4k_vasc2k_outputname6 = fullfile(dir_firstlevel,['w' subjectname '_' suffix '_atlas_2T4k_vasc2k_K1v.nii']);
                   twoT4k_vasc2k_outputname7 = fullfile(dir_firstlevel,['w' subjectname '_' suffix '_atlas_2T4k_vasc2k_k2v.nii']);
                   twoT4k_vasc2k_outputname8 = fullfile(dir_firstlevel,['w' subjectname '_' suffix '_atlas_2T4k_vasc2k_DV.nii']);
                   twoT4k_vasc2k_outputname9 = fullfile(dir_firstlevel,['w' subjectname '_' suffix '_atlas_2T4k_vasc2k_error.nii']);
                   LCN12_write_image(atlas_2T4k_vasc2k_K1,twoT4k_vasc2k_outputname1,'_' + suffix + '2T4k_vasc2k_K1',64,Vatlas_orig);        
                   LCN12_write_image(atlas_2T4k_vasc2k_k2,twoT4k_vasc2k_outputname2,'_' + suffix + '2T4k_vasc2k_k2',64,Vatlas_orig);        
                   LCN12_write_image(atlas_2T4k_vasc2k_k3,twoT4k_vasc2k_outputname3,'_' + suffix + '2T4k_vasc2k_k3',64,Vatlas_orig);        
                   LCN12_write_image(atlas_2T4k_vasc2k_k4,twoT4k_vasc2k_outputname4,'_' + suffix + '2T4k_vasc2k_k4',64,Vatlas_orig);        
                   LCN12_write_image(atlas_2T4k_vasc2k_Vb,twoT4k_vasc2k_outputname5,'_' + suffix + '2T4k_vasc2k_Vb',64,Vatlas_orig);  
                   LCN12_write_image(atlas_2T4k_vasc2k_K1v,twoT4k_vasc2k_outputname6,'_' + suffix + '2T4k_vasc2k_K1v',64,Vatlas_orig);
                   LCN12_write_image(atlas_2T4k_vasc2k_k2v,twoT4k_vasc2k_outputname7,'_' + suffix + '2T4k_vasc2k_k2v',64,Vatlas_orig);
                   LCN12_write_image(atlas_2T4k_vasc2k_DV,twoT4k_vasc2k_outputname8,'_' + suffix + '2T4k_vasc2k_DV',64,Vatlas_orig);        
                   LCN12_write_image(atlas_2T4k_vasc2k_error,twoT4k_vasc2k_outputname9,'_' + suffix + '2T4k_vasc2k_error',64,Vatlas_orig);
               end
               [AIC,SC] = LCN_calc_model_selection(nnz(WEIGHTS_PET),tmpB(index_min(1)).^2,nr_params_2T4k_Vb);
               results_region_based_AIC((subj-1)*4+4,2) = {modelname};
               results_region_based_SC((subj-1)*4+4,2)  = {modelname};
               results_region_based_AIC((subj-1)*4+4,2+voi) = {AIC};
               results_region_based_SC((subj-1)*4+4,2+voi)  = {SC};
               if figures_on == 1 || save_figures == 1 
                  figure(4)
                  subplot(224)
                  plot(MIDSCANTIMES,C_MEASURED,'ob:');
                  hold on
                  plot(MIDSCANTIMES(excluded_frames),C_MEASURED(excluded_frames),'*r'); % indicating unreliable data due to head movements
                  xlabel('time in min');
                  ylabel('kBq/ml');
                  title('2T4K\_vasc2k');
                  plot(MIDSCANTIMES,LCN_calc_model_2T4k_vasc2k(kfit_2T4k_vasc2k,MIDSCANTIMES,FRAMEDURATION,CA_FINETIMES,CA_WB_FINETIMES,FINETIMES,CALC_OPTION,STEP),'k')
                  axvalues = axis;
                  text(axvalues(2)*0.5,axvalues(3)+0.1*(axvalues(4)-axvalues(3)),['VD = ' num2str(kfit_2T4k_vasc2k(1)*(1+kfit_2T4k_vasc2k(3)/kfit_2T4k_vasc2k(4))/kfit_2T4k_vasc2k(2))]);
               end            
            end
            
        end % for loop models
        
        if save_figures == 1
           cd(dir_figures);
           saveas(hfig4,['fig_4models_' VOIname '_' suffix '.fig']);
        else 
           pause
        end   
                                
    end % for loop VOIs
    
    
    fprintf('atlas based analysis ended at %s \n',datetime('now'));
    fprintf(fid,'atlas based analysis ended at %s \n',datetime('now'));   
      

    % GENERATE A VOXEL-BASED PARAMETRIC IMAGE FOR LOGAN INPUT
    
    if do_voxelwise_logan
    
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

        % run voxel-based Logan analysis

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

        % save images

        outputname1 = fullfile(dir_firstlevel,['w' subjectname '_Logan_input_DV.nii']);
        outputname2 = fullfile(dir_firstlevel,['w' subjectname '_Logan_input_error.nii']);
        Vout1 = LCN12_write_image(DV_Logan,outputname1,'DV Logan',Vref.dt(1),Vref);        
        Vout2 = LCN12_write_image(DV_Logan_error,outputname2,'DV Logan error',Vref.dt(1),Vref); 

        fprintf('voxel based analysis for Logan_input ended at %s \n',datetime('now'));
        fprintf(fid,'voxel based analysis for Logan_input ended at %s \n',datetime('now'));
        
    end
    
end % for loop over subjects


%% SAVE THE RESULTS FOR THE REGION-BASED ANALYSIS
%--------------------------------------------------------------------------

for m = 1:nr_models
    modelname = model_list{m,1};
    if strcmp(modelname,'Logan_input')
       outputfile_excel = [outputfile_excel_region_based '_' suffix '_Logan_input.xlsx'];
       writetable(results_region_based_Logan_input_DV,outputfile_excel,'Sheet','DV');
       writetable(results_region_based_Logan_input_error,outputfile_excel,'Sheet','error');
    end
    if strcmp(modelname,'2T4k')
       outputfile_excel = [outputfile_excel_region_based '_' suffix '_2T4k.xlsx'];
       writetable(results_region_based_2T4k_K1,outputfile_excel,'Sheet','K1');
       writetable(results_region_based_2T4k_k2,outputfile_excel,'Sheet','k2');
       writetable(results_region_based_2T4k_k3,outputfile_excel,'Sheet','k3');
       writetable(results_region_based_2T4k_k4,outputfile_excel,'Sheet','k4');
       writetable(results_region_based_2T4k_DV,outputfile_excel,'Sheet','DV');
       writetable(results_region_based_2T4k_error,outputfile_excel,'Sheet','error');
    end
    if strcmp(modelname,'2T4k_Vb')
       outputfile_excel = [outputfile_excel_region_based '_' suffix '_2T4k_Vb.xlsx'];
       writetable(results_region_based_2T4k_Vb_K1,outputfile_excel,'Sheet','K1');
       writetable(results_region_based_2T4k_Vb_k2,outputfile_excel,'Sheet','k2');
       writetable(results_region_based_2T4k_Vb_k3,outputfile_excel,'Sheet','k3');
       writetable(results_region_based_2T4k_Vb_k4,outputfile_excel,'Sheet','k4');
       writetable(results_region_based_2T4k_Vb_Vb,outputfile_excel,'Sheet','Vb');
       writetable(results_region_based_2T4k_Vb_DV,outputfile_excel,'Sheet','DV');
       writetable(results_region_based_2T4k_Vb_error,outputfile_excel,'Sheet','error');
    end
    if strcmp(modelname,'2T4k_vasc1k')
       outputfile_excel = [outputfile_excel_region_based '_' suffix '_2T4k_vasc1k.xlsx'];
       writetable(results_region_based_2T4k_vasc1k_K1,outputfile_excel,'Sheet','K1');
       writetable(results_region_based_2T4k_vasc1k_k2,outputfile_excel,'Sheet','k2');
       writetable(results_region_based_2T4k_vasc1k_k3,outputfile_excel,'Sheet','k3');
       writetable(results_region_based_2T4k_vasc1k_k4,outputfile_excel,'Sheet','k4');
       writetable(results_region_based_2T4k_vasc1k_Vb,outputfile_excel,'Sheet','Vb');
       writetable(results_region_based_2T4k_vasc1k_K1v,outputfile_excel,'Sheet','K1v');
       writetable(results_region_based_2T4k_vasc1k_DV,outputfile_excel,'Sheet','DV');
       writetable(results_region_based_2T4k_vasc1k_error,outputfile_excel,'Sheet','error');
    end
    if strcmp(modelname,'2T4k_vasc2k')
       outputfile_excel = [outputfile_excel_region_based '_' suffix '_2T4k_vasc2k.xlsx'];
       writetable(results_region_based_2T4k_vasc2k_K1,outputfile_excel,'Sheet','K1');
       writetable(results_region_based_2T4k_vasc2k_k2,outputfile_excel,'Sheet','k2');
       writetable(results_region_based_2T4k_vasc2k_k3,outputfile_excel,'Sheet','k3');
       writetable(results_region_based_2T4k_vasc2k_k4,outputfile_excel,'Sheet','k4');
       writetable(results_region_based_2T4k_vasc2k_Vb,outputfile_excel,'Sheet','Vb');
       writetable(results_region_based_2T4k_vasc2k_K1v,outputfile_excel,'Sheet','K1v');
       writetable(results_region_based_2T4k_vasc2k_k2v,outputfile_excel,'Sheet','k2v');
       writetable(results_region_based_2T4k_vasc2k_DV,outputfile_excel,'Sheet','DV');
       writetable(results_region_based_2T4k_vasc2k_error,outputfile_excel,'Sheet','error');
    end
end

outputfile_excel = [outputfile_excel_region_based '_' suffix '_model_comparison_' char(datetime('today')) '.xlsx'];
writetable(results_region_based_AIC,outputfile_excel,'Sheet','AIC');
writetable(results_region_based_SC,outputfile_excel,'Sheet','SC');

outputfile_excel = [outputfile_excel_region_based '_' suffix '_nr_voxels_used_' char(datetime('today')) '.xlsx'];
writetable(results_number_voxels_VOI,outputfile_excel,'Sheet','nr of voxels');

cd(curdir);
fprintf('ALL DONE\n');