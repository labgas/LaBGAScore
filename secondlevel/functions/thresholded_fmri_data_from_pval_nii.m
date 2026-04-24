function [stat_img, fmri_dat, region_obj, region_table] = thresholded_fmri_data_from_pval_nii (fullpath_to_pval_nii, stat_vec, mask_obj, atlas_obj, p, p_type, stat_type, thr_type, k)

        stat_img = statistic_image(fullpath_to_pval_nii,'type','p');
        stat_img = apply_mask(stat_img,mask_obj);
        stat_img.sig = logical(stat_img.p < p); % manual thresholding at 0.05, statistic_image.threshold() breaks for some reason
        stat_img.p_type = p_type;
        
        if sum(stat_img.sig) == 0
            fprintf('\n No suprathreshold voxels after p threshold\n');
            fmri_dat = [];
            region_obj = [];
            region_table = [];
            
        else
        
            fmri_dat = fmri_data_st(stat_img);
            fmri_dat.dat = stat_vec;
            fmri_dat.dat_descrip = [stat_type ' thresholded at p_' thr_type ' < ' num2str(p)];
            path = fileparts(fullpath_to_pval_nii); 
            fmri_dat = fmri_dat.threshold([min(fmri_dat.dat(stat_img.sig))-0.001 max(fmri_dat.dat(stat_img.sig))+0.001],'raw-between','k',k);
            fmri_dat.image_names = [stat_type '_' strrep(num2str(p),'.','_') '_' thr_type '_k_' num2str(k) '.nii'];
            fmri_dat.fullpath = fullfile(path,[stat_type '_' strrep(num2str(p),'.','_') '_' thr_type '_k_' num2str(k) '.nii']);
            
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

