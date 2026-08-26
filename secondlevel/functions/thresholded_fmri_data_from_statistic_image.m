function [fmri_dat, region_obj, region_table] = thresholded_fmri_data_from_statistic_image (p_stat_img, stat_vec, atlas_obj, p, stat_type, thr_type, k)
% thresholded_fmri_data_from_statistic_image  Thresholded map and region table
% for a non-standard statistic, from an in-memory p-value statistic_image.
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
            fmri_dat.dat_descrip = [stat_type ' thresholded at p_' thr_type ' < ' num2str(p)];

            % Keep the SIGNIFICANT voxels, not every voxel whose statistic
            % happens to land in the significant range. Thresholding
            % 'raw-between' on [min(sig) max(sig)] retains anything inside that
            % window, which only equals masking when the p-value is monotonic
            % in the statistic. It is not: each voxel has its own permutation
            % null, so a non-significant voxel can easily carry a statistic
            % between the smallest and largest significant ones.
            % maskToSignificant parks the non-significant voxels on a sentinel
            % below the retained window, so the range threshold selects exactly
            % the significant set while 'k' still applies the extent rule.
            [dat_masked, lo, hi] = maskToSignificant(stat_vec, p_stat_img.sig);
            fmri_dat.dat = dat_masked;
            fmri_dat = fmri_dat.threshold([lo hi],'raw-between','k',k);
            fmri_dat.image_names = [stat_type '_' strrep(num2str(p),'.','_') '_' thr_type '_k_' num2str(k) '.nii'];
            
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

