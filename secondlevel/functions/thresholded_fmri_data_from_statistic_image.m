function [fmri_dat, region_obj, region_table] = thresholded_fmri_data_from_statistic_image (p_stat_img, stat_vec, atlas_obj, p, stat_type, thr_type, k)
%
% USAGE
%
% Create thresholded fmri_dat_st and region object for non-standard test statistic
% types (AUC, TFCE, etc) which can be used for visualization purposes
%
% 
% INPUTS
%
% - p_stat_img:     CANlab statistic_image object of type 'p' with p-values calculated on stat_vec thresholded at a value p of type thr_type
% - stat_vec:       vector of the length of p_stat_img.dat containing test statistic values on which the p-values in p_stat_img have been calculated
% - atlas_obj:      CANlab atlas object for labeling purposes
% - p:              scalar: threshold p-value used on p_stat_img
% - stat_type:      char array describing the type of test statistic in stat_vec
% - thr_type:       char array describing the type of correction done on p_stat_img, for example 'unc','fdr','fwe'
%                   only used for fmri_dat.fullpath and .dat_descrip here
% - k:              integer: additional extent-level threshold you may want to apply to fmri_dat, in voxels
%
% OUTPUTS
%
% -fmri_dat:        CANlab fmri_data_st object with stat_vec values in .dat
%                   thresholded at the threshold previously applied to the p_stat_img input (reflected in p and thr_type) 
%                   with additional extent-level threshold k
% -region_obj:      CANlab region object containing a region per significant cluster in fmri_dat for visualization purposes
% -region_table:    summary table summarizing the results in region_obj
%
% AUTHOR
% 
% Lukas Van Oudenhove, KU Leuven, 04/28/2023
        
        p_stat_img.sig = logical(p_stat_img.sig); % force .sig to logical

        if sum(p_stat_img.sig) == 0
            fprintf('\n No suprathreshold voxels after p threshold\n');
            fmri_dat = [];
            region_obj = [];
            region_table = [];
            
        else
        
            fmri_dat = fmri_data_st(p_stat_img);
            fmri_dat.dat = stat_vec;
            fmri_dat.dat_descrip = [stat_type ' thresholded at p_' thr_type ' < ' num2str(p)];
            fmri_dat = fmri_dat.threshold([min(fmri_dat.dat(p_stat_img.sig))-0.001 max(fmri_dat.dat(p_stat_img.sig))+0.001],'raw-between','k',k);
            fmri_dat.image_names = [stat_type '_' strrep(num2str(p),'.','_') '_' thr_type '_k_' num2str(k) '.nii'];
            
            if isempty(fmri_dat)
                fprintf('\n No suprathreshold voxels after k threshold\n');
                region_obj = [];
                region_table = [];
            
            else

                region_obj = region(fmri_dat);
                region_obj = region_obj.autolabel_regions_using_atlas(atlas_obj);

                [~,~,region_table] = table(region_obj);
                average_region = zeros(height(region_table),1);
                k_region = zeros(height(region_table),1);
                for r = 1:size(region_obj,2)
                    average_region(r,1) = mean(region_obj(1,r).val);
                    k_region(r,1) = region_obj(1,r).numVox;
                end
                region_table.minP = [];
                region_table.(['mean_' stat_type]) = average_region;
                region_table.k = k_region;
                
            end
            
        end
        
end

