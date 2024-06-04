%% LaBGAScore_atlas_binary_mask_from_atlas.m
%
%
% This script creates a binary mask by combining regions from one or more
% atlases, and automatically saves the fmri_mask_image, and .nii versions
% of it in maskdir of your dataset
%
% There are options to save the original atlas object created by
% select_atlas_subset (with one index for each individual parcel merged to
% create your mask), or an atlas object created from the fmri_mask_image
% object after merging the original parcels, i.e. with one index for the
% entire mask
% 
% USAGE
%
% Script should be run from the root directory of the superdataset, e.g.
% /data/proj_discoverie
%
% DEPENDENCIES
%
% CANlab's CanlabCore and Neuroimaging_Pattern_Masks Github repos on your Matlab path
% if needed, clone from https://github.com/canlab
%
% RESOURCES
%
% help atlas.select_atlas_subset
% https://canlab.github.io/_pages/using_canlab_atlases/using_canlab_atlases.html
%
% OPTIONS
%
% save_original_atlas_obj = true/false        saves original atlas object (i.e. BEFORE merging selected parcels into one fmri_mask_image object, hence one index per parcel) 
%                                               useful/needed if you want to label your parcel- or voxel-wise analyses with the labels included in your mask
%
% save_merged_atlas_obj = true/false          saves merged atlas object (i.e. AFTER merging selected parcels into one fmri_mask_image object, hence one index for the entire mask)
%                                               useful/needed if you want to extract roi averages 
%
% singleroi = true/false                      set to true if you are writing a mask/roi with one single contiguous region
%
%__________________________________________________________________________
%
% authors: Aleksandra Budzinska, Lukas Van Oudenhove
% date:   KU Leuven, July, 2022
%
%__________________________________________________________________________
% @(#)% LaBGAScore_atlas_binary_mask_from_atlas.m         v3.0       
% last modified: 2024/02/21


% Set options
% -------------------------------------------------------------------------

save_original_atlas_obj = true;
save_merged_atlas_obj = false;
single_roi = false;


% Define maskname and directory where mask will be written
%--------------------------------------------------------------------------

bit_rew_prep_s0_define_directories;
bit_rew_secondlevel_m1m_s0_a_set_up_paths_always_run_first;

modelname = 'bit_rew_m1m';

if single_roi

    roiname = 'reward_lateral_OFC_L';
    maskname = [modelname '_mask_' roiname];
    
else
    
    maskname_short = 'mask_reward_regions';
    maskname = [modelname '_' maskname_short];
    
end


% Load atlas of your choice
%--------------------------------------------------------------------------

atlasname1 = 'canlab2024_fine_2mm';
atlas1 = load_atlas(atlasname1); % you get a CANlab atlas object

%atlasname2 = 'cit168';
%atlas2 = load_atlas(atlasname2);

% ...

% NOTE: you can not only load CANlab atlases using keywords as in this example,
% but any atlas you like by creating a new atlas object from a .nii atlas image, 
% e.g. aicha = atlas('/opt/KUL_apps/spm12/atlas/AICHA.nii');
% see help atlas for adding labels etc


% Select regions you want to include from the lowest level of granularity
% ('labels_4')
%--------------------------------------------------------------------------

labels1 = {'insula_anterior','insula_operculum','Hippocampal_Formation','Amygdala','VStriatum','CAU','PUT','GP'};

atlas1_subset = select_atlas_subset(atlas1, labels1, 'labels_4');
atlas1_subset = atlas1_subset.threshold(0.20); % canlab2024 is a probabilistic atlas

atlas1_subset_flat = select_atlas_subset(atlas1, labels1, 'labels_4', 'flatten');
atlas1_subset_flat = atlas1_subset_flat.threshold(0.20);

% NOTE: you can select parcels using labels at different levels of
% granularity, in this case 'labels_4', which is a lot easier compared to
% having to list all the fine-grained labels!


% Select regions you want to include from the intermediate level of granularity
% ('labels_3')
%--------------------------------------------------------------------------

labels2 = {'vmPFC','Hythal','VTA_PBP','PAG','SN','LC+','Parabrachial'};

atlas2_subset = select_atlas_subset(atlas1, labels2, 'labels_3');
atlas2_subset = atlas2_subset.threshold(0.20); % canlab2024 is a probabilistic atlas

atlas2_subset_flat = select_atlas_subset(atlas1, labels2, 'labels_3', 'flatten');
atlas2_subset_flat = atlas2_subset_flat.threshold(0.20);


% Select regions you want to include from the highest level of granularity
% ('labels_2')
% -------------------------------------------------------------------------
labels3 = {'Ctx_47m_L','Ctx_a47r_L','Ctx_47s_L','Ctx_47m_R','Ctx_a47r_R','Ctx_47s_R','Ctx_a10p_L','Ctx_10pp_L','Ctx_11l_L','Ctx_13l_L','Ctx_OFC_L','Ctx_a10p_R','Ctx_10pp_R','Ctx_11l_R','Ctx_13l_R','Ctx_OFC_R'};

atlas3_subset = select_atlas_subset(atlas1, labels3, 'labels_2');
atlas3_subset = atlas3_subset.threshold(0.20); % canlab2024 is a probabilistic atlas

atlas3_subset_flat = select_atlas_subset(atlas1, labels3, 'labels_2','flatten');
atlas3_subset_flat = atlas3_subset_flat.threshold(0.20); % canlab2024 is a probabilistic atlas


% Change data type of probability maps if needed
% -------------------------------------------------------------------------

% NOTE: Before merging the atlases, make sure that they have the same type of variables (look at the probability_maps. e.g. single/sparse double)
% This command converts probability_maps property variable type from single to double, and then to sparse double - may not be needed for all atlases

%atlas2_subset.probability_maps = double(atlas2_subset.probability_maps);
%atlas2_subset.probability_maps = sparse(atlas2_subset.probability_maps);


% Merge the atlases
% -------------------------------------------------------------------------

atlases = {atlas1_subset, atlas2_subset, atlas3_subset};
atlases_flat = {atlas1_subset_flat, atlas2_subset_flat, atlas3_subset_flat};

atlas_nr = 1;

combined_atlas = atlases{1};
combined_atlas_flat = atlases_flat{1};
    
    while atlas_nr < nr_atlases
        
        combined_atlas = merge_atlases(combined_atlas, atlases{atlas_nr+1});
        combined_atlas_flat = merge_atlases(combined_atlas_flat, atlases_flat{atlas_nr+1});
        
        atlas_nr = atlas_nr+1;
        
    end
    
    if any(unique(combined_atlas_flat.dat) > 1)
        combined_atlas_flat.dat(combined_atlas_flat.dat > 0) = 1; % binarize
    end

% NOTE: merge_atlases function automatically resamples the space of the
% second to the first (which should not be needed if they are all subsets of the same atlas), and adds consecutive labels


% Convert the atlas object to a CANlab fmri_mask_image object
% -------------------------------------------------------------------------

mask = fmri_mask_image(combined_atlas_flat);
mask.volInfo_descrip = maskname;


% Create/save atlas objects according to options
% -------------------------------------------------------------------------

if save_original_atlas_obj
    
%     combined_atlas = atlas1_subset; % later scripts required a variable 'combined_atlas' in .mat file, comment this out if you have more than 1 subset
    
    if single_roi
        
        combined_atlas.atlas_name = roiname;
        
    else
        
        combined_atlas.atlas_name = maskname_short;
        
    end
    
end

if save_merged_atlas_obj
    
    roi_atlas = combined_atlas_flat;
    
    if single_roi
        
        roi_atlas.atlas_name = roiname;
        
    else
        
        roi_atlas.atlas_name = maskname_short;
        
    end
    
end

% Save your binarized mask in maskdir
%--------------------------------------------------------------------------

write(mask, 'fname', fullfile(maskdir,[maskname,'.nii'])); % writes to .nii file

savefilenamedata = fullfile(maskdir,[maskname,'.mat']); % saves fmri_mask_image object as .mat file

if save_original_atlas_obj && save_merged_atlas_obj
    
    save(savefilenamedata, 'combined_atlas','roi_atlas','mask', '-v7.3');
    fprintf('\nSaved mask & atlas objects\n');
    
    write(combined_atlas, 'fname', fullfile(maskdir,[maskname,'_atlas.nii'])); % writes mask to .nii file
    
elseif save_original_atlas_obj && ~save_merged_atlas_obj
    
    save(savefilenamedata, 'combined_atlas','mask', '-v7.3');
    fprintf('\nSaved mask & atlas object\n');
    
    write(combined_atlas, 'fname', fullfile(maskdir,[maskname,'_atlas.nii'])); % writes mask to .nii file
    
elseif ~save_original_atlas_obj && save_merged_atlas_obj
    
    save(savefilenamedata, 'roi_atlas','mask', '-v7.3');
    fprintf('\nSaved mask & atlas object\n');
    
else 
    
    save(savefilenamedata, 'mask', '-v7.3');
    fprintf('\nSaved mask\n');
    
end
