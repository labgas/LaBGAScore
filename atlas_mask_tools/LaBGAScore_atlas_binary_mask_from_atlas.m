%% LaBGAScore_atlas_binary_mask_from_atlas_ROI.m
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
% @(#)% LaBGAScore_atlas_binary_mask_from_atlas.m         v2.0       
% last modified: 2024/02/15


% Set options
% -------------------------------------------------------------------------

save_original_atlas_obj = false;
save_merged_atlas_obj = true;
single_roi = true;


% Define maskname and directory where mask will be written
%--------------------------------------------------------------------------

ery_4b_prep_s0_define_directories;
ery_4b_secondlevel_m1m_s0_a_set_up_paths_always_run_first;

modelname = 'ery_4b_m1m';

if single_roi

    roiname = 'reward_lateral_OFC_L';
    maskname = [modelname '_mask_' roiname];
    
else
    
    maskname_short = 'mask_reward_all_regions';
    maskname = [modelname '_' maskname_short];
    
end


% Load atlas of your choice
%--------------------------------------------------------------------------

atlasname1 = 'canlab2018';
atlas1 = load_atlas(atlasname1); % you get a CANlab atlas object

%atlasname2 = 'cit168';
%atlas2 = load_atlas(atlasname2);

% ...

% NOTE: you can not only load CANlab atlases using keywords as in this example,
% but any atlas you like by creating a new atlas object from a .nii atlas image, 
% e.g. aicha = atlas('/opt/KUL_apps/spm12/atlas/AICHA.nii');
% see help atlas for adding labels etc


% Select regions you want to include from the first atlas 
%--------------------------------------------------------------------------

labels1 = {'Ctx_47m_L','Ctx_a47r_L','Ctx_47s_L'};
atlas1_subset = select_atlas_subset(atlas1, labels1);

% NOTE: you can also use the numbers rather than the labels of regions to
% select as a vector
% see help atlas.select_atlas_subset

% NOTE: the code commented out below is for the case where you want to
% combine parcels from different atlases, which is deprecated since we will
% use the excellent canlab2023 atlas from now on


% Select regions you want to include from the second atlas
% -------------------------------------------------------------------------
%labels2 = {};
%atlas2_subset = select_atlas_subset(atlas2, labels2); % still a CANlab atlas object


% Change data type of probability maps if needed
% -------------------------------------------------------------------------

% NOTE: Before merging the atlases, make sure that they have the same type of variables (look at the probability_maps. e.g. single/sparse double)
% This command converts probability_maps property variable type from single to double, and then to sparse double - may not be needed for all atlases

%atlas2_subset.probability_maps = double(atlas2_subset.probability_maps);
%atlas2_subset.probability_maps = sparse(atlas2_subset.probability_maps);


% Merge the atlases
% -------------------------------------------------------------------------

%combined_atlas = merge_atlases(atlas1_subset, atlas2_subset); 

% NOTE: merge_atlases function automatically resamples the space of the second to the first, and adds consecutive labels


% Convert the atlas object to a CANlab fmri_mask_image object
% -------------------------------------------------------------------------

mask = fmri_mask_image(atlas1_subset);
mask.volInfo_descrip = maskname;

% NOTE: this automatically binarizes the mask!
% TO BE CHECKED IF THIS IS ALSO THE CASE FOR PROBABILISTIC ATLASES


% Create/save atlas objects according to options
% -------------------------------------------------------------------------

if save_original_atlas_obj
    
    combined_atlas = atlas1_subset; % later scripts required a variable 'combined_atlas' in .mat file, comment this out if you have more than 1 subset
    
    if single_roi
        
        combined_atlas.atlas_name = roiname;
        
    else
        
        combined_atlas.atlas_name = maskname_short;
        
    end
    
end

if save_merged_atlas_obj
    
    if single_roi 
        
        roi_atlas = atlas(mask);
        
        roi_atlas.atlas_name = roiname;
        roi_atlas.labels = {labels1};
        roi_atlas.label_descriptions = {roiname};
        
    else
        
        fprintf('\n');
        warning('option to save merged atlas object not yet implemented for multiple non-contiguous rois');
        fprintf('\n');
        
% NOTE: the code below may work but needs to be tested thoroughly before it
% can possibly be implemented!
% Alternative would be to define each rois as a separate atlas subset in
% the code above, convert them all to fmri_mask_image objects, then to
% region objects - this could all be looped. 
% Finally, those region objects can be concatenated and converted into a
% single atlas object using region2atlas().
%
%         roi_region = region(mask);
%         roi_atlas = 
%         roi_atlas.atlas_name = maskname_short;
%         roi_atlas.label_descriptions = cell(1,1); % to add appropriate labels and their descriptions here in case multiple non-contiguous rois are included in the mask/atlas, each roi would need to be defined as a separate subset above
        
    end
    
end

% Save your binarized mask in maskdir
%--------------------------------------------------------------------------

write(mask, 'fname', fullfile(maskdir,[maskname,'.nii'])); % writes to .nii file

savefilenamedata = fullfile(maskdir,[maskname,'.mat']); % saves fmri_mask_image object as .mat file

if save_original_atlas_obj && save_merged_atlas_obj
    
    save(savefilenamedata, 'combined_atlas','roi_atlas','mask', '-v7.3');
    fprintf('\nSaved mask & atlas objects\n');
    
elseif save_original_atlas_obj && ~save_merged_atlas_obj
    
    save(savefilenamedata, 'combined_atlas','mask', '-v7.3');
    fprintf('\nSaved mask & atlas object\n');
    
elseif ~save_original_atlas_obj && save_merged_atlas_obj
    
    save(savefilenamedata, 'roi_atlas','mask', '-v7.3');
    fprintf('\nSaved mask & atlas object\n');
    
else 
    
    save(savefilenamedata, 'mask', '-v7.3');
    fprintf('\nSaved mask\n');
    
end
