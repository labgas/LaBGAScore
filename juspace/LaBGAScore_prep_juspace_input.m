%% LaBGAScore_prep_juspace_input.m
%
%
% *USAGE*
%
% This script creates the necessary input files (i.e. firstlevel con images) 
% for spatial correlation analysis with PET receptor maps with the JuSpace toolbox 
% from LaBGAS standard first and second level file organization
%
% Launch JuSpace from your Matlab terminal by typing JuSpace
%
% For more info about JuSpace, see the following resources
%
% # JuSpace Github repo: https://github.com/juryxy/JuSpace
%
% # JuSpace HBM paper: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7814756/
%
%
% *OPTIONS*
%
% No options need to be specified in this script, but it uses GLM options
% set in the a2_set_default_options.m CANlab_help_examples (LaBGAS fork)
% script for the corresponding second-level model, particularly
%
% *atlasname_glm
% *maskname_glm
%
% Help prep_3a_run_second_level_regression_and_save for more info about
% these options
%
%
% -------------------------------------------------------------------------
%
% author: Lukas Van Oudenhove
%
% date:   KU Leuven, August, 2023
%
% -------------------------------------------------------------------------
%
% LaBGAScore_prep_juspace_input.m         v2.0
%
% last modified: 2023/11/21
%
%
%% PREP WORK
% -------------------------------------------------------------------------

% GET MODEL-SPECIFIC PATHS AND OPTIONS

a_set_up_paths_always_run_first;


% LOAD DSGN, DAT, AND DATA OBJECTS

if ~exist('DSGN','var') || ~exist('DAT','var')
    
    load(fullfile(resultsdir,'image_names_and_setup.mat'));
    
end

if ~exist('DATA_OBJ','var') || ~exist('DATA_OBJsc','var')
    
    load(fullfile(resultsdir,'data_objects.mat'));
    load(fullfile(resultsdir,'data_objects_scaled.mat'));
    
end


% GET FIRSTLEVEL DIR INFO

firstmodeldir = DSGN.modeldir;

firstlist = dir(fullfile(firstmodeldir,'sub-*'));
firstsubjs = cellstr(char(firstlist(:).name));
firstsubjdirs = cell(size(firstsubjs,1),1);

    for firstsub = 1:size(firstsubjs,1)
        firstsubjdirs{firstsub,1} = fullfile(firstlist(firstsub).folder,firstlist(firstsub).name);
    end
    
clear firstlist


% SET JUSPACE DIRS AND CREATE IF NEEDED

juspacedir = fullfile(resultsdir,'juspace');
if ~isfolder(juspacedir)
    mkdir(juspacedir);
end

juspaceinputdir = fullfile(juspacedir, 'input_files');
if ~isfolder(juspaceinputdir)
    mkdir(juspaceinputdir);
end

juspaceatlasdir = fullfile(juspaceinputdir, 'atlas');
if ~isfolder(juspaceatlasdir)
    mkdir(juspaceatlasdir);
end

juspacemaskdir = fullfile(juspaceinputdir, 'mask');
if ~isfolder(juspacemaskdir)
    mkdir(juspacemaskdir);
end

juspaceresultsdir = fullfile(juspacedir, 'results');
if ~isfolder(juspaceresultsdir)
    mkdir(juspaceresultsdir);
end

for cond = 1:size(DAT.conditions,2)
    if ~isfolder(fullfile(juspaceresultsdir,DAT.conditions{cond}))
        mkdir(fullfile(juspaceresultsdir,DAT.conditions{cond}));
    end
end

for cont = 1:size(DAT.contrastnames,2)
    if ~isfolder(fullfile(juspaceresultsdir,DAT.contrastnames{cont}))
        mkdir(fullfile(juspaceresultsdir,DAT.contrastnames{cont}));
    end
end

for sub = 1:size(firstsubjs,1)
    if ~isfolder(fullfile(juspaceinputdir,firstsubjs{sub}))
        mkdir(fullfile(juspaceinputdir,firstsubjs{sub}));
    end
end


% ATLAS PREP WORK

if exist('atlasname_glm','var') && ~isempty(atlasname_glm) && exist(atlasname_glm, 'file')

    if contains(atlasname_glm,'.mat')
        
        load(atlasname_glm);
        [~, atlasname_short] = fileparts(atlasname_glm);
        
        combined_atlas.atlas_name = atlasname_short;
        combined_atlas.image_names = [atlasname_short '.mat'];
        combined_atlas.fullpath = atlasname_glm;
        
        write(combined_atlas,'fname',fullfile(juspaceatlasdir,[atlasname_short '.nii']),'overwrite');
        
        fprintf('\nUse custom-made atlas %s in JuSpace\n\n', atlasname_short);

    elseif ischar(atlasname_glm)

        combined_atlas = load_atlas(atlasname_glm);
        
        if contains(atlasname_glm,'canlab2023')
                combined_atlas = combined_atlas.threshold(0.25); % only keep probability values > 0.25
        end
        
        combined_atlas.atlas_name = atlasname_glm;
        combined_atlas.image_names = [atlasname_glm '.nii'];
        combined_atlas.fullpath = which([atlasname_glm '.nii']);
        
        write(combined_atlas,'fname',fullfile(juspaceatlasdir,[atlasname_glm '.nii']),'overwrite');
        
        fprintf('\nUse custom-made atlas %s in JuSpace\n\n', atlasname_glm);

    else

         error('\ninvalid option "%s" defined in atlasname_glm variable, should be a keyword for load_atlas.m or a .mat file containing an atlas object, check docs\n\n',atlasname_glm)

    end
    
else
    
    fprintf('\nUse standard Neuromorphometrics atlas provided with JuSpace, or any other atlas file of your choice in NIfTI format \n\n');
    
end


% MASK PREP WORK

clear mask

if exist('maskname_glm','var') && ~isempty(maskname_glm) && exist(maskname_glm, 'file')

    mask_img = fmri_mask_image(maskname_glm,'noverbose');
    target = fmri_data(fullfile(firstsubjdirs{1},'con_0001.nii'),'noverbose'); % resample to space of functional images
    voxelsize_target = abs(diag(target.volInfo.mat(1:3, 1:3)))';
    voxelsize_glmmask = abs(diag(mask_img.volInfo.mat(1:3, 1:3)))';
    
    if ~isequal(voxelsize_glmmask,voxelsize_cat_obj)
        mask2write = resample_space(mask_img,target);
        mask2write.dat(mask2write.dat > 0) = 1; % binarize mask
    else
        mask2write = mask_img;
        mask2write.dat(mask2write.dat > 0) = 1; % binarize mask
    end
    
    maskname_final = [maskname_glm(1:end-4) '_binary'];
    [~, maskname_short] = fileparts(maskname_final);
    mask = fullfile(juspacemaskdir,[maskname_short '.nii']);
    write(mask2write,'fname',mask,'overwrite');
    
    fprintf('\nMasking con images with custom mask %s \n\n', maskname_short);
    
else
    
    default_mask = which('gm_mask_canlab2023_coarse_fmriprep20_0_25.nii');
    mask_img = fmri_mask_image(default_mask,'noverbose');
    target = fmri_data(fullfile(firstsubjdirs{1},'con_0001.nii'),'noverbose'); % resample to space of functional images
    voxelsize_target = abs(diag(target.volInfo.mat(1:3, 1:3)))';
    voxelsize_default_mask = abs(diag(mask_img.volInfo.mat(1:3, 1:3)))';
    
    if ~isequal(voxelsize_glmmask,voxelsize_cat_obj)
        mask2write = resample_space(mask_img,target);
        mask2write.dat(mask2write.dat > 0) = 1; % binarize mask
    else
        mask2write = mask_img;
        mask2write.dat(mask2write.dat > 0) = 1; % binarize mask
    end
    
    maskname_final = [default_mask(1:end-4) '_binary'];
    [~, maskname_short] = fileparts(maskname_final);
    mask = fullfile(juspacemaskdir,[maskname_short '.nii']);
    write(mask2write,'fname',mask,'overwrite');
    
    fprintf('\nMasking con images with default mask %s \n\n', maskname_short);
    
end


%% MASK CON IMAGES AND MOVE TO JUSPACE INPUT DIRS
%--------------------------------------------------------------------------

kc = size(DAT.conditions, 2);

for c = 1:kc
    
    fprintf('\n\n');
    printhdr(['CONDITION #', num2str(c), ': ', upper(DAT.conditions{c})]);
    fprintf('\n\n');
    
    % SELECT DATA FOR THIS CONDITION
    
    switch myscaling_glm

                case 'raw'
                    fprintf('\nRaw (unscaled) condition images used for spatial correlations with PET receptor maps\n\n');
                    scaling_string = 'no_scaling';
                    cat_obj = DATA_OBJ{c};

                case 'scaled'
                    fprintf('\nZ-scored condition images used for spatial correlations with PET receptor maps\n\n');
                    scaling_string = 'scaling_z_score_conditions';
                    cat_obj = DATA_OBJsc{c};

                otherwise
                    error('\nInvalid option "%s" defined in myscaling_glm variable in a2_set_default_options script, choose between "raw",  and "scaled"\n\n', myscaling_glm);

    end
    
    % MASK CON IMAGES 
    
    cat_obj_masked = apply_mask(cat_obj,mask2write);
    
    % WRITE MASKED CON IMAGES
    
    cat_obj_2write = cat_obj_masked;
    cat_obj_2write.fullpath = string(cat_obj_2write.fullpath);
    firstmodelparts = strsplit(firstmodeldir,'/');
    newstr = [firstmodelparts{end-1} '/' firstmodelparts{end}];
    juspaceparts = strsplit(juspaceinputdir,'/');
    newstr2 = [juspaceparts{end-4} '/' juspaceparts{end-3} '/' juspaceparts{end-2} '/' juspaceparts{end-1} '/' juspaceparts{end}];
    cat_obj_2write.fullpath = replace(cat_obj_2write.fullpath,newstr,newstr2);
    cat_obj_2write.fullpath = replace(cat_obj_2write.fullpath,['con_000' num2str(c)],DAT.conditions{c});
    cat_obj_2write.fullpath = char(cat_obj_2write.fullpath);
    
    for i = 1:size(cat_obj_2write.fullpath,1)
        obj_2write = get_wh_image(cat_obj_2write,i);
        write(obj_2write)
    end
    
end % for loop over contrasts
