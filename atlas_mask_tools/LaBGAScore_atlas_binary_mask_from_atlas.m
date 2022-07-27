%% LaBGAScore_atlas_binary_mask_from_atlas.m
%
% This script creates a binary mask by combining regions from one or more
% atlases, and automatically saves the atlas, fmri_mask_image, and .nii versions
% of it in maskdir of your dataset
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
%__________________________________________________________________________
%
% authors: Aleksandra Budzinska, Lukas Van Oudenhove
% date:   KU Leuven, July, 2022
%
%__________________________________________________________________________
% @(#)% LaBGAScore_atlas_binary_mask_from_atlas.m         v1.1       
% last modified: 2022/07/26


% Load the atlas
canlab = load_atlas('canlab2018'); % you get a CANlab atlas object

% Select regions you want to include from the selected atlas 
canlab_subset = select_atlas_subset(canlab, {'Ctx_p24pr_L','Ctx_p24pr_R','Ctx_33pr_L','Ctx_33pr_R','Ctx_a24pr_L','Ctx_a24pr_R','Ctx_p32pr_L','Ctx_p32pr_R','Ctx_a24_L','Ctx_a24_R','Ctx_d32_L','Ctx_d32_R','Ctx_8BM_L','Ctx_8BM_R','Ctx_p32_L','Ctx_p32_R','Ctx_10r_L','Ctx_10r_R','Ctx_47m_L','Ctx_47m_R','Ctx_10d_L','Ctx_10d_R','Ctx_a47r_L','Ctx_a47r_R','Ctx_a10p_L','Ctx_a10p_R','Ctx_10pp_L','Ctx_10pp_R','Ctx_11l_L','Ctx_11l_R','Ctx_13l_L','Ctx_13l_R','Ctx_OFC_L','Ctx_OFC_R','Ctx_47s_L','Ctx_47s_R','Ctx_52_L','Ctx_52_R','Ctx_PFcm_L','Ctx_PFcm_R','Ctx_PoI2_L','Ctx_PoI2_R','Ctx_FOP4_L','Ctx_FOP4_R','Ctx_MI_L','Ctx_MI_R','Ctx_Pir_L','Ctx_Pir_R','Ctx_AVI_L','Ctx_AVI_R','Ctx_AAIC_L','Ctx_AAIC_R','Ctx_FOP1_L','Ctx_FOP1_R','Ctx_FOP3_L','Ctx_FOP3_R','Ctx_FOP2_L','Ctx_FOP2_R','Ctx_25_L','Ctx_25_R','Ctx_s32_L','Ctx_s32_R','Ctx_pOFC_L','Ctx_pOFC_R','Ctx_PoI1_L','Ctx_PoI1_R','Ctx_Ig_L','Ctx_Ig_R','Ctx_FOP5_L','Ctx_FOP5_R','Ctx_p10p_L','Ctx_p10p_R','Ctx_PI_L','Ctx_PI_R','Ctx_a32pr_L','Ctx_a32pr_R','Ctx_p24_L','Ctx_p24_R','Thal_VPM','Bstem_PAG','Bstem_SC_R','Bstem_SC_L','Bstem_IC_R','IC_L','Dorsal_raphe_DR','Median_raphe_MR_R','pbn_R','pbn_L','lc_R','lc_L','rvm_R','nrm','dmnx_nts_R','dmnx_nts_L','ncf_R','ncf_L','ncs_B6_B8','nrp_B5','nuc_ambiguus_R','medullary_raphe','spinal_trigeminal_R','spinal_trigeminal_L','CA2_Hippocampus_','Amygdala_CM_','DG_Hippocampus_','Amygdala_SF_','CA3_Hippocampus_','Amygdala_AStr_','Amygdala_LB_','CA1_Hippocampus_','Subiculum'});

% Load and selecet the regions from another atlas (e.g. subcortical) like above 
subcortical = load_atlas('cit168'); % you get a CANlab atlas object
subcortical_subset = select_atlas_subset(subcortical, {'Put','Cau','NAC','BST_SLEA','GPe','GPi','SNc','RN','SNr','PBP','VTA','VeP','Hythal','STN'}); % still a CANlab atlas object

% Before merging the atlases, make sure that they have the same type of variables (look at the probability_maps. E.g. single/sparse double)
% This command converts probability_maps property variable type from single to double, and then to sparse double - may not be needed for all atlases
subcortical_subset.probability_maps = double(subcortical_subset.probability_maps);
subcortical_subset.probability_maps = sparse(subcortical_subset.probability_maps);

% Merge the atlases
combined_atlas = merge_atlases(subcortical_subset, canlab_subset); % merge_atlases function automatically resamples the space of the second to the first, and adds consecutive labels

% Binarize the merged atlases by converting the atlas object to a CANlab fmri_mask_image object 
mask = fmri_mask_image(combined_atlas);

% Save your binarized mask as a NIFTI file
write(mask, 'fname', 'C:\Users\lukas\Downloads\mask.nii', 'overwrite');