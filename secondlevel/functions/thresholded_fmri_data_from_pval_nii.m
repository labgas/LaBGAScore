function [stat_img, fmri_dat, region_obj, region_table] = thresholded_fmri_data_from_pval_nii (fullpath_to_pval_nii, stat_vec, mask_obj, atlas_obj, p, p_type, stat_type, thr_type, k)
% thresholded_fmri_data_from_pval_nii  Thresholded map and region table for a
% non-standard statistic, from a p-value NIfTI on disk.
%
% USAGE
%
% Create thresholded fmri_dat_st and region object for non-standard test statistic
% types (AUC, TFCE, etc) which can be used for visualization purposes
%
% INPUTS
%
% - fullpath_to_pval_nii:       fullpath to a .nii file with p-values calculated on stat_vec thresholded at a value p of type thr_type
% - stat_vec:                   vector containing test statistic values on which the p-values in fullpath_to_pval_nii have been calculated
% - mask_obj:                   CANlab mask_image object for masking fullpath_to_pval_nii
% - atlas_obj:                  CANlab atlas object for labeling purposes
% - p:                          scalar: threshold p-value used on fullpath_to_pval_nii
% - p_type:                     char array describing the method used to calculated p-values in fullpath_to_pval_nii, for example 'based on 1000 permutations' 
% - stat_type:                  char array describing the type of test statistic in stat_vec
% - thr_type:                   char array describing the type of correction done on p_stat_img, for example 'unc','fdr','fwe'
%                               only used for fmri_dat.fullpath and .dat_descrip here
% - k:                          integer: additional extent-level threshold you may want to apply to fmri_dat, in voxels
%
% OUTPUTS
%
% -fmri_dat:                    CANlab fmri_data_st object with stat_vec values in .dat
%                               thresholded at the threshold previously applied to the fullpath_to_pval_nii input (reflected in p, p_type, and thr_type) 
%                               with additional extent-level threshold k
% -region_obj:                  CANlab region object containing a region per significant cluster in fmri_dat for visualization purposes
% -region_table:                summary table summarizing the results in region_obj
%
% AUTHOR
% 
% Lukas Van Oudenhove, KU Leuven, 04/28/2023

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
            fmri_dat.dat_descrip = [stat_type ' thresholded at p_' thr_type ' < ' num2str(p)];
            out_dir = fileparts(fullpath_to_pval_nii);

            % See maskToSignificant: thresholding 'raw-between' on the range of
            % the significant values also retains non-significant voxels whose
            % statistic falls inside that window, which happens whenever the
            % p-value is not monotonic in the statistic (it is not, because each
            % voxel has its own permutation null).
            [dat_masked, lo, hi] = maskToSignificant(stat_vec, stat_img.sig);
            fmri_dat.dat = dat_masked;
            fmri_dat = fmri_dat.threshold([lo hi],'raw-between','k',k);
            fmri_dat.image_names = [stat_type '_' strrep(num2str(p),'.','_') '_' thr_type '_k_' num2str(k) '.nii'];
            fmri_dat.fullpath = fullfile(out_dir,[stat_type '_' strrep(num2str(p),'.','_') '_' thr_type '_k_' num2str(k) '.nii']);
            
            if isempty(fmri_dat)
                fprintf('\n No suprathreshold voxels after k threshold\n');
                region_obj = [];
                region_table = [];
            
            else

                region_obj = region(fmri_dat);
                region_obj = region_obj.autolabel_regions_using_atlas(atlas_obj);

                [~,~,region_table] = table(region_obj);

                % Preallocate and loop over the SAME count. These used to
                % disagree: preallocation used height(region_table) while the
                % loop ran over size(region_obj,2), so if those differed the
                % array grew silently and the write-back below failed with a
                % table-height error.
                nReg = numel(region_obj);
                average_region = zeros(nReg,1);
                k_region = zeros(nReg,1);
                for r = 1:nReg
                    average_region(r,1) = mean(region_obj(r).val);
                    k_region(r,1) = region_obj(r).numVox;
                end

                if height(region_table) ~= nReg
                    % Returning the plain table beats dying in a visualisation
                    % helper, but say so clearly.
                    warning('thresholded_fmri_data:regionCountMismatch', ...
                        ['region() returned %d regions but its table has %d rows, so the ' ...
                         'mean-%s and k columns cannot be appended. Returning the ' ...
                         'un-augmented table.'], nReg, height(region_table), stat_type);
                else
                    % Assigning [] to a table variable DELETES it, and deleting
                    % a column that does not exist is an error rather than a
                    % no-op, so check first.
                    if ismember('minP', region_table.Properties.VariableNames)
                        region_table.minP = [];
                    end
                    region_table.(['mean_' stat_type]) = average_region;
                    region_table.k = k_region;
                end
                
            end
            
        end
        
end

