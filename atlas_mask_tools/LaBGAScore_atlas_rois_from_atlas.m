%% LaBGAScore_atlas_rois_from_atlas.m
%
%
% This script creates a rois by combining regions from one or more atlases, 
% and automatically saves the rois as atlas objects, and .nii versions
% of it in maskdir of your dataset
%
% There are options to save the original atlas object created by
% select_atlas_subset (with one index for each individual parcel merged to
% create your rois)
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
% save_original_roi_atlas_obj = true/false      saves original roi atlas objects (i.e. BEFORE merging selected parcels into one roi atlas object, hence one index per parcel) 
%                                               useful/needed if you want to label your voxel-wise analyses with the labels included in your mask
%
%                                               flat roi atlas objects (i.e. AFTER merging selected parcels into one roi atlas object, hence one index for the entire roi)
%                                               will always be saved as they are needed to extract roi averages using functionality in prep_3a script
%
% write_roi2nii = true/false                    writes .nii files for each of the flat roi atlas objects, which can be useful for some purposes
%                                               
%
%__________________________________________________________________________
%
% authors: Aleksandra Budzinska, Lukas Van Oudenhove
% date:   KU Leuven, July, 2022
%
%__________________________________________________________________________
% @(#)% LaBGAScore_atlas_rois_from_atlas.m         v1.4       
% last modified: 2025/05/07


%% SET OPTIONS, DEFINE DIRS, AND SPECIFY MODEL AND ROI SET NAMES
% -------------------------------------------------------------------------

save_original_roi_atlas_obj = true;
write_roi2nii = true;

% STUDY- AND MODEL-SPECIFIC SCRIPTS/VARS

LaBGAScore_prep_s0_define_directories;
a_set_up_paths_always_run_first;

roi_modelname = 'bit_rew_m1m';          % prefix which will be added to names of saved roi files which will be written in model-specific maskdir
roi_set_name = 'reward_regions';        % descriptive name for set of rois which will be included in filename


%% LOAD ATLAS(ES)
%--------------------------------------------------------------------------

atlasname1 = 'canlab2024_fine_2mm';
atlas1 = load_atlas(atlasname1); % CANlab atlas object

%atlasname2 = 'cit168';
%atlas2 = load_atlas(atlasname2);

% ...

% NOTE: you can not only load CANlab atlases using keywords as in this example,
% but any atlas you like by creating a new atlas object from a .nii atlas image, 
% e.g. aicha = atlas('/opt/KUL_apps/spm12/atlas/AICHA.nii');
% see help atlas for adding labels etc


%% SELECT REGIONS TO BE INCLUDED AT THE LOWEST LEVEL OF GRANULARITY
%--------------------------------------------------------------------------

% NOTE: this is level 'labels_4' in canlab2024

granularity1 = 'labels_4';
labels1 = {{'insula_anterior','insula_operculum'},{'VStriatum'},{'CAU'},{'PUT'}}; % cell array of cell arrays grouping labels that should be joined into one ROI
roi_names1 = {'amINS','ventral_striatum','caudate','putamen'};
roi_bilateral_idx1 = [true true true true];

nr_rois1 = ((size(roi_names1,2)*2) - (size(roi_names1,2) - sum(roi_bilateral_idx1)));
roi_atlases1_flat = cell(1, nr_rois1);
roi_counter = 1;

if save_original_roi_atlas_obj
    
    roi_atlases1 = cell(1, nr_rois1);
    
end

while roi_counter < nr_rois1

    for label = 1:size(labels1,2)

            if roi_bilateral_idx1(label)

                roi_atlas_bilateral = select_atlas_subset(atlas1,labels1{label}, granularity1);

                roi_atlases1_flat{roi_counter} = select_atlas_subset(roi_atlas_bilateral.threshold(0.20),{'_L'},'flatten');
                roi_atlases1_flat{roi_counter}.atlas_name = [roi_names1{label} '_L'];
                
                if save_original_roi_atlas_obj
                    
                    roi_atlases1{roi_counter} = select_atlas_subset(roi_atlas_bilateral.threshold(0.20),{'_L'});
                    roi_atlases1{roi_counter}.atlas_name = [roi_names1{label} '_L'];
                    
                end

                    roi_counter = roi_counter + 1;

                roi_atlases1_flat{roi_counter} = select_atlas_subset(roi_atlas_bilateral.threshold(0.20),{'_R'},'flatten');
                roi_atlases1_flat{roi_counter}.atlas_name = [roi_names1{label} '_R'];
                
                if save_original_roi_atlas_obj
                    
                    roi_atlases1{roi_counter} = select_atlas_subset(roi_atlas_bilateral.threshold(0.20),{'_R'});
                    roi_atlases1{roi_counter}.atlas_name = [roi_names1{label} '_R'];
                    
                end

                    roi_counter = roi_counter + 1;

            else

                roi_atlases1_flat{roi_counter} = select_atlas_subset(atlas1.threshold(0.20),labels1{label}, granularity1, 'flatten');
                roi_atlases1_flat{roi_counter}.atlas_name = roi_names1{label};
                
                if save_original_roi_atlas_obj
                    
                    roi_atlases1{roi_counter} = select_atlas_subset(atlas1.threshold(0.20),labels1{label}, granularity1);
                    roi_atlases1{roi_counter}.atlas_name = roi_names1{label};
                    
                end
                    

                    roi_counter = roi_counter + 1;

            end
            
            clear roi_atlas_bilateral
            
    end

end

clear roi_counter


%% SELECT REGIONS TO BE INCLUDED AT THE INTERMEDIATE LEVEL OF GRANULARITY
%--------------------------------------------------------------------------

% NOTE: this is level 'labels_3' in canlab2024

granularity2 = 'labels_3';
labels2 = {{'vmPFC'},{'Hythal'},{'VTA_PBP'}}; % cell array of cell arrays grouping labels that should be joined into one ROI
roi_names2 = {'vmPFC','hypothalamus','VTA'};
roi_bilateral_idx2 = [true false false];

nr_rois2 = ((size(roi_names2,2)*2) - (size(roi_names2,2) - sum(roi_bilateral_idx2)));
roi_atlases2_flat = cell(1, nr_rois2);
roi_counter = 1;

if save_original_roi_atlas_obj
    
    roi_atlases2 = cell(1, nr_rois2);
    
end

while roi_counter < nr_rois2

    for label = 1:size(labels2,2)

            if roi_bilateral_idx2(label)

                roi_atlas_bilateral = select_atlas_subset(atlas1,labels2{label}, granularity2);

                roi_atlases2_flat{roi_counter} = select_atlas_subset(roi_atlas_bilateral.threshold(0.20),{'_L'},'flatten');
                roi_atlases2_flat{roi_counter}.atlas_name = [roi_names2{label} '_L'];
                
                if save_original_roi_atlas_obj
                    
                    roi_atlases2{roi_counter} = select_atlas_subset(roi_atlas_bilateral.threshold(0.20),{'_L'});
                    roi_atlases2{roi_counter}.atlas_name = [roi_names2{label} '_L'];
                    
                end

                    roi_counter = roi_counter + 1;

                roi_atlases2_flat{roi_counter} = select_atlas_subset(roi_atlas_bilateral.threshold(0.20),{'_R'},'flatten');
                roi_atlases2_flat{roi_counter}.atlas_name = [roi_names2{label} '_R'];
                
                if save_original_roi_atlas_obj
                    
                    roi_atlases2{roi_counter} = select_atlas_subset(roi_atlas_bilateral.threshold(0.20),{'_R'});
                    roi_atlases2{roi_counter}.atlas_name = [roi_names2{label} '_R'];
                    
                end

                    roi_counter = roi_counter + 1;

            else

                roi_atlases2_flat{roi_counter} = select_atlas_subset(atlas1.threshold(0.20),labels2{label}, granularity2, 'flatten');
                roi_atlases2_flat{roi_counter}.atlas_name = roi_names2{label};
                
                if save_original_roi_atlas_obj
                    
                    roi_atlases2{roi_counter} = select_atlas_subset(atlas1.threshold(0.20),labels2{label}, granularity2);
                    roi_atlases2{roi_counter}.atlas_name = roi_names2{label};
                    
                end

                    roi_counter = roi_counter + 1;

            end
            
            clear roi_atlas_bilateral
            
    end

end

clear roi_counter


%% SELECT REGIONS TO BE INCLUDED AT THE HIGHEST LEVEL OF GRANULARITY
%--------------------------------------------------------------------------

% NOTE: this is level 'labels_2' in canlab2024

granularity3 = 'labels_2';
labels3 = {{'Ctx_47m_L','Ctx_a47r_L','Ctx_47s_L','Ctx_47m_R','Ctx_a47r_R','Ctx_47s_R'},{'Ctx_a10p_L','Ctx_10pp_L','Ctx_11l_L','Ctx_13l_L','Ctx_OFC_L','Ctx_a10p_R','Ctx_10pp_R','Ctx_11l_R','Ctx_13l_R','Ctx_OFC_R'}}; % cell array of cell arrays grouping labels that should be joined into one ROI
roi_names3 = {'lOFC','mOFC'};
roi_bilateral_idx3 = [true true];

nr_rois3 = ((size(roi_names3,2)*2) - (size(roi_names3,2) - sum(roi_bilateral_idx3)));
roi_atlases3_flat = cell(1, nr_rois3);
roi_counter = 1;

if save_original_roi_atlas_obj
    
    roi_atlases3 = cell(1, nr_rois3);
    
end

while roi_counter < nr_rois3

    for label = 1:size(labels3,2)

            if roi_bilateral_idx3(label)

                roi_atlas_bilateral = select_atlas_subset(atlas1,labels3{label}, granularity3);

                roi_atlases3_flat{roi_counter} = select_atlas_subset(roi_atlas_bilateral.threshold(0.20),{'_L'},'flatten');
                roi_atlases3_flat{roi_counter}.atlas_name = [roi_names3{label} '_L'];
                
                if save_original_roi_atlas_obj
                    
                    roi_atlases3{roi_counter} = select_atlas_subset(roi_atlas_bilateral.threshold(0.20),{'_L'});
                    roi_atlases3{roi_counter}.atlas_name = [roi_names3{label} '_L'];
                    
                end

                    roi_counter = roi_counter + 1;

                roi_atlases3_flat{roi_counter} = select_atlas_subset(roi_atlas_bilateral.threshold(0.20),{'_R'},'flatten');
                roi_atlases3_flat{roi_counter}.atlas_name = [roi_names3{label} '_R'];
                
                if save_original_roi_atlas_obj
                    
                    roi_atlases3{roi_counter} = select_atlas_subset(roi_atlas_bilateral.threshold(0.20),{'_R'});
                    roi_atlases3{roi_counter}.atlas_name = [roi_names3{label} '_R'];
                    
                end

                    roi_counter = roi_counter + 1;

            else

                roi_atlases3_flat{roi_counter} = select_atlas_subset(atlas1.threshold(0.20),labels3{label}, granularity3, 'flatten');
                roi_atlases3_flat{roi_counter}.atlas_name = roi_names3{label};
                
                if save_original_roi_atlas_obj
                    
                    roi_atlases3{roi_counter} = select_atlas_subset(atlas1.threshold(0.20),labels3{label}, granularity3);
                    roi_atlases3{roi_counter}.atlas_name = roi_names3{label};
                    
                end

                    roi_counter = roi_counter + 1;

            end
            
            clear roi_atlas_bilateral
            
    end

end

clear roi_counter


% Merge the cell arrays with atlas objects
% -------------------------------------------------------------------------

roi_atlases_flat = [roi_atlases1_flat,roi_atlases2_flat,roi_atlases3_flat];

    if save_original_roi_atlas_obj
        
       roi_atlases = [roi_atlases1,roi_atlases2,roi_atlases3];
       
    end


% Save your roi atlas objects in maskdir
%--------------------------------------------------------------------------

savefilenamedata = fullfile(maskdir,[roi_modelname '_rois_' roi_set_name '.mat']); % saves cell array with roi atlas objects as .mat file

if save_original_roi_atlas_obj
    
    save(savefilenamedata, 'roi_atlases', 'roi_atlases_flat', '-v7.3');
    fprintf('\nSaved roi atlas objects\n');
    
else 
    
    save(savefilenamedata, 'roi_atlases_flat', '-v7.3');
    fprintf('\nSaved roi atlas objects\n');
    
end

if write_roi2nii
    
    for r = 1:size(roi_atlases_flat,2)
        
        saveroinamedata = fullfile(maskdir,[roi_atlases_flat{r}.atlas_name '.nii']);
        write(roi_atlases_flat{r},'fname',saveroinamedata);
        clear saveroinamedata
        
    end
    
end


% Create montage of ROIs
%--------------------------------------------------------------------------

roi_atlas = roi_atlases_flat{1};
roi_atlas.atlas_name = 'combined_roi_atlas';
roi_atlas.labels{1} = roi_atlases_flat{1}.atlas_name;
      
roi = 2;
    
    while roi < (size(roi_atlases_flat,2) + 1)
        
        roi_atlas = merge_atlases(roi_atlas, roi_atlases_flat{roi},'noreplace');
%           LVO: whether or not 'noreplace' option is chosen or not should
%           not matter if all regions are from same atlas

        roi_atlas.labels{roi} = roi_atlases_flat{roi}.atlas_name;
        
        roi = roi+1;
        
    end
    
cmap2cell = scn_standard_colors(size(roi_atlas.labels,2));
cmap2 = zeros(size(roi_atlas.labels,2),3);

    for color = 1:size(cmap2cell,2)
        cmap2(color,:) = cmap2cell{color};
    end

figure;
o2 = canlab_results_fmridisplay([], 'compact');
o2 = addblobs(o2, atlas2region(roi_atlas),'indexmap',cmap2,'interp','nearest');
o2 = title_montage(o2, 5, 'atlas used for extraction of roi averages');
set(gcf,'WindowState','maximized');
drawnow,snapnow;

clear o2

