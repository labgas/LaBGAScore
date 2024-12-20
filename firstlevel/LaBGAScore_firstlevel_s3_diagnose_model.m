%% LaBGAScore_firstlevel_s3_diagnose_model.m
%
% This script will run diagnostics on first level models and publish them
% as an html report using the CANlab function scn_spm_design_check, and
% save variance inflation factors as a .mat file
%
%
% *DOCUMENTATION*
% 
% type "help scn_spm_design_check" in Matlab command window
%
%
% *USAGE*
%
% Script should be called from LaBGAScore_firstlevel_s2_fit_model.m, it
% is not for standalone use
%
% -------------------------------------------------------------------------
%
% authors: Lukas Van Oudenhove
%
% date:   March, 2022
%
% -------------------------------------------------------------------------
%
% LaBGAScore_firstlevel_s3_diagnose_model.m         v1.4   
%
% last modified: 2024/08/07
%
%
%% VIFS AND DESIGN
% -------------------------------------------------------------------------

    vifs = scn_spm_design_check(subjfirstdir,'events_only','vif_thresh',LaBGAS_options.mandatory.vif_thresh);
    drawnow, snapnow
    save('vifs','vifs');
    close all
    

%% FIRST LEVEL MAPS
% -------------------------------------------------------------------------

if LaBGAS_options.display.plotmontages

    load(fullfile(subjfirstdir,'SPM.mat'));

    tmapnames = dir(fullfile(subjfirstdir,'spmT_*.nii'));
    
        if isfield(DSGN,'singletrials')
            if ~isempty(DSGN.singletrials)
                tmapnames = tmapnames(1:size(DSGN.contrastnames,2));
            end
        end
    
    tmapspaths = cell(1,size(tmapnames,1));
    tmapsobj = cell(1,size(tmapnames,1));
    montages = cell(1,size(tmapnames,1));

    [~,maskname,maskext] = fileparts(LaBGAS_options.display.mask);
    mask = [maskname,maskext];
    mask_img = fmri_mask_image(LaBGAS_options.display.mask,'noverbose');
    target = fullfile(tmapnames(1).folder,tmapnames(1).name);
    mask_img = resample_space(mask_img,target); % resample mask space to T-map space
    mask_img.dat(mask_img.dat > 0) = 1; % binarize resampled mask

    fprintf('\nShowing results at p < %s %s, k = %s, mask = %s\n',num2str(LaBGAS_options.display.input_threshold),LaBGAS_options.display.thresh_type,num2str(LaBGAS_options.display.k),mask); 

        for tmap = 1:size(tmapnames,1)
            tmapspaths{tmap} = fullfile(tmapnames(tmap).folder,tmapnames(tmap).name);
            tmapsobj{tmap} = statistic_image('image_names',tmapspaths{tmap},'type','T','dfe',SPM.xX.erdf);
            tmapsobj{tmap} = threshold(tmapsobj{tmap},LaBGAS_options.display.input_threshold,LaBGAS_options.display.thresh_type,'k',LaBGAS_options.display.k,'mask',mask_img);
            create_figure('fmridisplay');
            wh_montage = 5;
            axis off
            figtitle = DSGN.contrastnames{tmap};
            o3 = canlab_results_fmridisplay([],'outline','linewidth',0.5,'montagetype','compact','overlay','mni_icbm152_t1_tal_nlin_sym_09a_brainonly.img');
            o3 = addblobs(o3,tmapsobj{tmap},'splitcolor',{[.1 .8 .8] [.1 .1 .8] [.9 .4 0] [1 1 0]});
            [o3,title_handle] = title_montage(o3,wh_montage,figtitle);
            set(title_handle,'FontSize',18);
            fighan = activate_figures(o3);
            f3 = fighan{1};
            f3.Tag = figtitle;
            f3.WindowState = 'maximized';
            drawnow,snapnow
            close(f3)
            clear figtitle o3 title_handle fighan f3
        end
        
end