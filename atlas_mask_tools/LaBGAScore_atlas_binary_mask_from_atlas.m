%% LaBGAScore_atlas_binary_mask_from_atlas.m
%
%
% *USAGE*
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
% Script should be run from the root directory of the superdataset, e.g.
% /data/proj_discoverie
%
%
% *OPTIONS*
%
% * save_original_atlas_obj = true/false        saves original atlas object (i.e. BEFORE merging selected parcels into one fmri_mask_image object, hence one index per parcel)
%                                                 useful/needed if you want to label your parcel- or voxel-wise analyses with the labels included in your mask
%
% * save_merged_atlas_obj = true/false          saves merged atlas object (i.e. AFTER merging selected parcels into one fmri_mask_image object, hence one index for the entire mask)
%                                                 useful/needed if you want to extract roi averages
%
% * single_roi = true/false                     set to true if you are writing a mask/roi with one single contiguous region
%
% * prob_threshold = <0-1>                      probability threshold applied to each selected subset; canlab2024 and
%                                                 other probabilistic atlases need this. Was hardcoded as 0.20 in six
%                                                 separate places before v3.3
%
% * strict_labels = true/false                  true (default) errors when a requested label matches nothing in the
%                                                 atlas; false downgrades that to a warning and builds the mask from
%                                                 whatever did match. Leaving a level EMPTY ({}) is always legitimate
%                                                 and is skipped silently - the check is for typos, which would
%                                                 otherwise produce a quietly smaller mask
%
%
% *DEPENDENCIES*
%
% CANlab's CanlabCore and Neuroimaging_Pattern_Masks Github repos on your Matlab path
% if needed, clone from https://github.com/canlab
%
%
% *NOTES*
%
% help atlas.select_atlas_subset
% https://canlab.github.io/_pages/using_canlab_atlases/using_canlab_atlases.html
%
% GRANULARITY LEVELS MAY BE LEFT EMPTY. Set any of labels1/labels2/labels3 to
% {} and that level is skipped; the remaining levels are merged as usual. This
% works for the FIRST level too, which the pre-v3.3 code could not do because
% it seeded the merge from atlases{1} unconditionally. Before v3.3 an empty
% level died at .threshold() with 'Empty or invalid .probability_maps. Cannot
% threshold', because select_atlas_subset returns a 0-region atlas when its
% loop over the requested labels never executes.
%
% This script is checked in as a worked example for one specific study/model
% (reward-region atlas labels, maskname 'bit_rew_m1m', etc.) rather than a fully
% generic template - the atlas/label/ROI names in the body below are
% study-specific and need adapting for your own study even though the
% directory-setup calls above have been genericized
%
% -------------------------------------------------------------------------
%
% modified by: Aleksandra Budzinska, Lukas Van Oudenhove
%
% date:   KU Leuven, July, 2022
%
% -------------------------------------------------------------------------
%
% LaBGAScore_atlas_binary_mask_from_atlas.m         v3.3
%
% last modified: 2026/11/26


% Set options
% -------------------------------------------------------------------------

save_original_atlas_obj = true;
save_merged_atlas_obj = false;
single_roi = false;
prob_threshold = 0.20;
strict_labels = true;


% Define maskname and directory where mask will be written
%--------------------------------------------------------------------------

LaBGAScore_prep_s0_define_directories;
a_set_up_paths_always_run_first;

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


% Select the regions you want to include, per level of granularity
%--------------------------------------------------------------------------

% NOTE: you can select parcels using labels at different levels of
% granularity rather than having to list every fine-grained label - in
% canlab2024 'labels_4' is the coarsest and 'labels_2' the finest.
%
% Leave a level EMPTY ({}) if your model does not draw on it. Empty levels
% are skipped and the remaining ones merged as usual - including when it is
% the FIRST level that is empty. A level that is NOT empty but matches
% nothing in the atlas is an error rather than a silently smaller mask; see
% strict_labels in the options above.

% lowest level of granularity
labels1 = {'insula_anterior','insula_operculum','Hippocampal_Formation','Amygdala','VStriatum','CAU','PUT','GP'};
granularity1 = 'labels_4';

% intermediate level of granularity
labels2 = {'vmPFC','Hythal','VTA_PBP','PAG','SN','LC+','Parabrachial'};
granularity2 = 'labels_3';

% highest level of granularity
labels3 = {'Ctx_47m_L','Ctx_a47r_L','Ctx_47s_L','Ctx_47m_R','Ctx_a47r_R','Ctx_47s_R','Ctx_a10p_L','Ctx_10pp_L','Ctx_11l_L','Ctx_13l_L','Ctx_OFC_L','Ctx_a10p_R','Ctx_10pp_R','Ctx_11l_R','Ctx_13l_R','Ctx_OFC_R'};
granularity3 = 'labels_2';


% Build one thresholded subset per NON-EMPTY level
% -------------------------------------------------------------------------

level_labels = {labels1, labels2, labels3};
level_granularity = {granularity1, granularity2, granularity3};

atlases = {};
atlases_flat = {};
atlas_fields = fieldnames(atlas1);

for level = 1:numel(level_labels)

    these_labels = level_labels{level};
    this_granularity = level_granularity{level};

    if isempty(these_labels)

        fprintf('\nLevel %d (%s): no labels requested - skipping this level\n', level, this_granularity);
        continue

    end

    if ~any(strcmp(this_granularity, atlas_fields))

        error('Level %d: ''%s'' is not a field of atlas ''%s''. Available label fields: %s', ...
            level, this_granularity, atlasname1, ...
            strjoin(atlas_fields(contains(atlas_fields, 'labels'))', ', '));

    end

    % select_atlas_subset matches SUBSTRINGS unless 'exact'/'regexp' is
    % passed, so test each requested label the same way rather than guessing
    available_labels = atlas1.(this_granularity);
    is_unmatched = cellfun(@(s) ~any(contains(available_labels, s)), these_labels);

    if any(is_unmatched)

        msg = sprintf(['Level %d (%s): %d of %d requested labels match nothing in atlas ''%s'':\n  %s\n' ...
            'Check the spelling and the granularity level. Set strict_labels = false to downgrade this to a warning.'], ...
            level, this_granularity, sum(is_unmatched), numel(these_labels), atlasname1, ...
            strjoin(these_labels(is_unmatched), ', '));

        if strict_labels
            error('%s', msg);
        else
            warning('%s', msg);
        end

    end

    subset = select_atlas_subset(atlas1, these_labels, this_granularity);

    if num_regions(subset) == 0

        error(['Level %d (%s): the requested labels selected 0 atlas regions, so .threshold() ' ...
            'would fail with ''Empty or invalid .probability_maps''. Labels requested: %s'], ...
            level, this_granularity, strjoin(these_labels, ', '));

    end

    subset = subset.threshold(prob_threshold); % canlab2024 is a probabilistic atlas

    subset_flat = select_atlas_subset(atlas1, these_labels, this_granularity, 'flatten');
    subset_flat = subset_flat.threshold(prob_threshold);

    atlases{end+1} = subset; %#ok<SAGROW>
    atlases_flat{end+1} = subset_flat; %#ok<SAGROW>

    fprintf('\nLevel %d (%s): %d labels requested, %d atlas regions selected\n', ...
        level, this_granularity, numel(these_labels), num_regions(subset));

end

if isempty(atlases)

    error(['All granularity levels are empty, so there is nothing to build a mask from. ' ...
        'Populate at least one of labels1/labels2/labels3.']);

end


% Change data type of probability maps if needed
% -------------------------------------------------------------------------

% NOTE: Before merging the atlases, make sure that they have the same type of variables (look at the probability_maps. e.g. single/sparse double)
% This command converts probability_maps property variable type from single to double, and then to sparse double - may not be needed for all atlases

%atlases{2}.probability_maps = double(atlases{2}.probability_maps);
%atlases{2}.probability_maps = sparse(atlases{2}.probability_maps);


% Merge the atlases
% -------------------------------------------------------------------------

combined_atlas = atlases{1};
combined_atlas_flat = atlases_flat{1};

for atlas_nr = 2:numel(atlases)

    combined_atlas = merge_atlases(combined_atlas, atlases{atlas_nr});
    combined_atlas_flat = merge_atlases(combined_atlas_flat, atlases_flat{atlas_nr});

end

fprintf('\nMerged %d granularity level(s) into a combined atlas with %d regions\n', ...
    numel(atlases), num_regions(combined_atlas));

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
