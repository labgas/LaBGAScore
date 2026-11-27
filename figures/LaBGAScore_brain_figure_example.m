%% LaBGAScore_brain_figure_example.m
%
%
% *USAGE*
%
% Worked examples of LaBGAScore_brain_figure, one section per layout, that
% double as its test: run the whole script and every layout is drawn and
% exported. Each section is self-contained, so copy the one you need into
% your study's code and swap in your own map.
%
% 1. load CanlabCore's bundled emotionreg sample (Wager et al. 2008, n = 30)
%   and derive a t-map, a thresholded t-map, a Bayes-factor map, an atlas
%   subset and a mask from it
%   NOTE: the sample is used so the script runs anywhere CanlabCore is
%       installed. Section 16 shows where the same objects live in the
%       results files of the LaBGAS second-level templates
%
% 2. one section per layout: the CanlabCore preset montages, custom slice
%   rows, unthresholded maps with significance outlines, multi-threshold
%   layers, Bayes factors, accurately projected cortical surfaces, cutaways,
%   subcortical scenes, per-region close-ups, zoom insets (magnified, or as
%   topographies of boxes, clusters or atlas ROIs), atlases, small
%   multiples, the two topography modes, and a multi-panel composite
%
% Every figure is written to outdir as png at dpi; add 'pdf' / 'svg' to
% 'formats' for the versions a journal asks for.
%
% Run it headless, which is the LaBGAS default:
%   matlab -nodisplay -r "run('LaBGAScore_brain_figure_example.m'); exit" < /dev/null
%
%
% *OPTIONS*
%
% * outdir      where the figures are written (default: a folder under
%               tempdir). NOT the current directory: run() changes into the
%               script's own folder while it runs, so a pwd-based default would
%               write into the LaBGAScore repo - or, from a study's code
%               subdataset, into a DataLad dataset
% * dpi         resolution of the png files (default 300 here, to keep the
%               gallery quick; the function's own default is 600)
%
%
% *DEPENDENCIES*
%
% CanlabCore and Neuroimaging_Pattern_Masks (atlases) on the Matlab path
% SPM12 (image i/o)
% Curve Fitting Toolbox, for the 'roi' topography only (tpaps)
%
%
% *NOTES*
%
% emotionreg is in MNI152NLin6Asym space (SPM normalisation), so the surface
% sections pass 'sourcespace', 'MNI152NLin6Asym'. LaBGAS data are fmriprep 20
% output in MNI152NLin2009cAsym, the function's default: leave 'sourcespace'
% out for your own maps. Its 3.4 x 3.4 x 4.5 mm voxels, in a different
% template from the default underlay, also make its blobs overhang the
% cortical edge more than LaBGAS 2 mm data will; see 'clip_mask'.
%
% See figures/README.md for the design decisions behind the defaults.
%
% -------------------------------------------------------------------------
%
% author: Lukas Van Oudenhove & Claude Opus 5.5
%
% date:   KU Leuven, November, 2026
%
% -------------------------------------------------------------------------
%
% LaBGAScore_brain_figure_example.m         v1.0
%
% last modified: 2026/11/26
%
%
%% SETTINGS
%--------------------------------------------------------------------------

outdir = fullfile(tempdir, 'LaBGAScore_brain_figure_examples');
dpi = 300;

if ~exist(outdir, 'dir'), mkdir(outdir); end
fprintf('Writing the example figures to %s\n', outdir);
sample_space = {'sourcespace', 'MNI152NLin6Asym'};   % emotionreg only; see NOTES


%% 1. LOAD THE SAMPLE DATA AND DERIVE THE MAPS
%--------------------------------------------------------------------------

dat = load_image_set('emotionreg', 'noverbose');         % fmri_data, 30 contrast images
t = ttest(dat);                                          % statistic_image, unthresholded
t_fdr = threshold(t, .05, 'fdr', 'noverbose');           % thresholded copy
BF = estimateBayesFactor(t, 't');                        % statistic_image of type 'BF', stored as 2*ln(BF10)

atl = load_atlas('canlab2024_coarse_fmriprep20_2mm');
reward_atl = select_atlas_subset(atl, {'VStriatum', 'Amygdala', 'insula'}, 'labels_3');
reward_mask = fmri_mask_image(reward_atl);               % one binary mask


%% 2. PRESET MONTAGE (canlab_results_fmridisplay 'compact2')
%--------------------------------------------------------------------------
% CanlabCore's own layout, with this function's colours, size and export.

fig = LaBGAScore_brain_figure(t_fdr, 'layout', 'compact2', ...
    'save', fullfile(outdir, '02_compact2'), 'dpi', dpi);
close(fig);


%% 3. CUSTOM SLICE ROWS, SLICE LOCATOR AND LABELS
%--------------------------------------------------------------------------
% One row mixing a sagittal and axials; the locator marks where the axials
% are, and the strongest regions are named from the atlas.

fig = LaBGAScore_brain_figure(t_fdr, 'layout', 'slices', ...
    'slices', {{{'sagittal', -2}, {'axial', -20:10:50}}}, ...
    'slice_locator', true, 'labels', 'auto', 'max_labels', 6, ...
    'save', fullfile(outdir, '03_slices_labels'), 'dpi', dpi);
close(fig);


%% 4. UNTHRESHOLDED MAP WITH SIGNIFICANCE OUTLINES (Acil et al. 2026 Fig 2)
%--------------------------------------------------------------------------
% The whole map pale; FDR-significant voxels at full colour, outlined. For a
% map without p-values (e.g. bootstrapped SVM weights) pass the significant
% voxels as 'sig' instead: a thresholded statistic_image, region or mask.

fig = LaBGAScore_brain_figure(t, 'layout', 'slices', ...
    'style', 'unthresholded_outline', 'sig', {.05, 'fdr'}, ...
    'slices', {{'axial', -20:10:50}}, ...
    'save', fullfile(outdir, '04_unthresholded_outline'), 'dpi', dpi);
close(fig);


%% 5. MULTI-THRESHOLD LAYERS (Spisak et al. 2026 Fig 3)
%--------------------------------------------------------------------------
% Lenient level first; each level gets its own colours and legend entry.

fig = LaBGAScore_brain_figure(t, 'layout', 'slices', ...
    'thresholds', {{.01, 'unc'}, {.05, 'fdr'}}, ...
    'slices', {{'axial', -20:10:50}}, ...
    'save', fullfile(outdir, '05_multi_threshold'), 'dpi', dpi);
close(fig);


%% 6. BAYES FACTORS (Kwon et al. 2026 Fig 2a)
%--------------------------------------------------------------------------
% Red = evidence for H1 (BF10 >= 3), blue = evidence for H0 (BF10 <= 1/3),
% colorbar in BF10. The map is stored as 2*ln(BF10); info.bayes reports the
% range in BF10 units, and a note when strong evidence for H0 is out of reach
% at this sample size.

[fig, ~, info] = LaBGAScore_brain_figure(BF, 'layout', 'slices', ...
    'slices', {{'axial', -20:10:50}}, ...
    'save', fullfile(outdir, '06_bayes'), 'dpi', dpi);
disp(info.bayes)
close(fig);


%% 7. CORTICAL SURFACES, ACCURATELY PROJECTED
%--------------------------------------------------------------------------
% fsLR_32k inflated surfaces with MNIsurf projection (Wu et al. 2018).

fig = LaBGAScore_brain_figure(t_fdr, 'layout', 'surfaces', sample_space{:}, ...
    'save', fullfile(outdir, '07_surfaces'), 'dpi', dpi);
close(fig);


%% 8. CUTAWAY: INSULA SLAB (Kwon et al. 2026 Fig 1b)
%--------------------------------------------------------------------------

fig = LaBGAScore_brain_figure(t_fdr, 'layout', 'cutaway', 'cutaway', 'left_insula_slab', ...
    'width', 'single', 'save', fullfile(outdir, '08_insula_slab'), 'dpi', dpi);
close(fig);


%% 9. SUBCORTICAL SCENE
%--------------------------------------------------------------------------

fig = LaBGAScore_brain_figure(t_fdr, 'layout', 'subcortex', ...
    'structures', {'thalamus', 'caudate', 'putamen', 'amygdala', 'brainstem'}, ...
    'save', fullfile(outdir, '09_subcortex'), 'dpi', dpi);
close(fig);


%% 10. ONE CLOSE-UP PER REGION
%--------------------------------------------------------------------------
% At most 'max_regions' (20 by default, the LaBGAS report limit), each named
% from the atlas at its peak voxel.

fig = LaBGAScore_brain_figure(t_fdr, 'layout', 'regioncenters', ...
    'threshold', {.05, 'fdr', 'k', 20}, 'labels', 'auto', 'max_regions', 12, ...
    'save', fullfile(outdir, '10_regioncenters'), 'dpi', dpi);
close(fig);


%% 11. ZOOM INSETS (Kragel et al. 2023 Fig 1a)
%--------------------------------------------------------------------------
% Boxes are [u1 u2 v1 v2] in mm in the slice plane: (y, z) for a sagittal.

fig = LaBGAScore_brain_figure(t_fdr, 'layout', 'zoom', 'width', 'single', ...
    'zoom', {{'sagittal', 38, [-5 35 -25 15; -45 -5 -5 35], {'anterior insula', 'posterior insula'}}}, ...
    'save', fullfile(outdir, '11_zoom'), 'dpi', dpi);
close(fig);


%% 11b. TOPOGRAPHY INSETS ZOOMING OUT OF CLUSTERS (Kragel et al. 2023 Fig 1)
%--------------------------------------------------------------------------
% The 'max_zoom' strongest clusters, each on an axial slice through its peak,
% boxed, and above it the UNTHRESHOLDED values of the boxed plane as a 3-D
% surface. All insets share one scale and their own colorbar.

fig = LaBGAScore_brain_figure(t_fdr, 'layout', 'zoom', 'zoom', 'clusters', ...
    'zoom_inset', 'topography', ...
    'save', fullfile(outdir, '11b_zoom_topography_clusters'), 'dpi', dpi);
close(fig);


%% 11c. TOPOGRAPHY INSETS ZOOMING OUT OF ATLAS ROIs
%--------------------------------------------------------------------------
% A coronal slice through each ROI's centre, the ROI outlined; the inset is
% the spline-fitted ROI topography over the slice's own axes (x, z), so both
% hemispheres show side by side.

fig = LaBGAScore_brain_figure(t_fdr, 'layout', 'zoom', 'zoom', {'VStriatum', 'Amygdala'}, ...
    'orientation', 'coronal', 'zoom_inset', 'topography', ...
    'save', fullfile(outdir, '11c_zoom_topography_rois'), 'dpi', dpi);
close(fig);


%% 11d. TOPOGRAPHY INSETS FROM YOUR OWN BOXES
%--------------------------------------------------------------------------
% Section 11 with 'zoom_inset', 'topography'.

fig = LaBGAScore_brain_figure(t_fdr, 'layout', 'zoom', 'width', 'single', ...
    'zoom', {{'sagittal', 38, [-5 35 -25 15; -45 -5 -5 35], {'anterior insula', 'posterior insula'}}}, ...
    'zoom_inset', 'topography', ...
    'save', fullfile(outdir, '11d_zoom_topography_boxes'), 'dpi', dpi);
close(fig);


%% 12. AN ATLAS IN UNIQUE COLOURS, WITH A LEGEND
%--------------------------------------------------------------------------
% The subset holds 36 fine parcels - too many for a legend, which is only
% drawn up to 24. 'atlas_level' colours them by their labels_3 name instead:
% 10 entries, VStriatum / Amygdala / three insula zones, left and right.

fig = LaBGAScore_brain_figure(reward_atl, 'layout', 'slices', 'atlas_level', 'labels_3', ...
    'slices', {{'axial', -20:8:20}}, ...
    'save', fullfile(outdir, '12_atlas'), 'dpi', dpi);
close(fig);


%% 13. A MASK OUTLINE OVER A MAP
%--------------------------------------------------------------------------

fig = LaBGAScore_brain_figure(t_fdr, 'layout', 'slices', ...
    'overlay', {{reward_mask, 'outline', [0 0 0]}}, ...
    'slices', {{'axial', -20:8:20}}, ...
    'save', fullfile(outdir, '13_mask_outline'), 'dpi', dpi);
close(fig);


%% 14. SMALL MULTIPLES: ONE PANEL PER PARTICIPANT (Sun et al. 2026 Fig 3)
%--------------------------------------------------------------------------
% Shared colour limits across panels; a ring marks the same coordinate in each.
% 'image' picks which images: pass the whole object and a list, rather than
% subsetting it with get_wh_image, which on an fmri_data_st (DATA_OBJ_CON in
% the LaBGAS templates) returned all images instead of the ones asked for.

fig = LaBGAScore_brain_figure(dat, 'layout', 'small_multiples', 'image', 1:9, ...
    'panel', {'axial', 0}, 'ncols', 3, 'width', 'single', ...
    'titles', arrayfun(@(k) sprintf('Participant %d', k), 1:9, 'UniformOutput', false), ...
    'markers', [-30 20 0], ...
    'save', fullfile(outdir, '14_small_multiples'), 'dpi', dpi);
close(fig);


%% 15. TOPOGRAPHY: CLUSTER SURFACE PLOTS (Kragel et al. 2023 Fig 1c)
%--------------------------------------------------------------------------
% 'roi': voxel values over |x| and y within the ventral striatum, fitted with
% a thin-plate smoothing spline. info.panels{1}.topography holds the
% smoothing parameter tpaps chose and descriptive slopes per mm.
% 'slice': the raw values of one plane in a box - the general form of
% figures/cluster_surface_plots.m, in mm instead of voxel indices.

[fig, ~, info] = LaBGAScore_brain_figure(t, 'layout', 'topography', 'roi', {'VStriatum'}, ...
    'zlabel', 't', 'width', 'single', 'save', fullfile(outdir, '15a_topography_roi'), 'dpi', dpi);
disp(info.panels{1}.topography)
close(fig);

fig = LaBGAScore_brain_figure(t, 'layout', 'topography', 'topo_mode', 'slice', ...
    'plane', {'sagittal', 38}, 'box', [-10 30 -20 20], ...
    'zlabel', 't', 'width', 'single', 'save', fullfile(outdir, '15b_topography_slice'), 'dpi', dpi);
close(fig);


%% 16. A MULTI-PANEL PAPER FIGURE
%--------------------------------------------------------------------------
% Slices, surfaces and a subcortical scene share one colour mapping and one
% colorbar; a topography panel of the same contrast keeps its own. Panel
% letters are added automatically. Export pdf as well for the journal.

fig = LaBGAScore_brain_figure(t_fdr, 'layout', 'composite', sample_space{:}, ...
    'panels', { ...
        {'slices', 'slices', {{'axial', -20:10:50}}, 'labels', 'auto', 'max_labels', 5}, ...
        {'surfaces'}, ...
        {'subcortex'}, ...
        {'topography', 'map', t, 'roi', {'VStriatum'}, 'zlabel', 't'}}, ...
    'grid', [3 2], 'save', fullfile(outdir, '16_composite'), 'formats', {'png', 'pdf'}, 'dpi', dpi);
close(fig);


%% 17. WHERE THESE OBJECTS LIVE IN LaBGAS RESULTS FILES
%--------------------------------------------------------------------------
% The second-level templates of the CANlab_help_examples LaBGAS fork save the
% maps this function takes. With c the contrast number and resultsdir the
% model's results directory:
%
%   load(fullfile(resultsdir, 'regression_stats_and_maps_<...>.mat'), 'regression_stats_results');
%   t  = regression_stats_results{c}.t;     % statistic_image, unthresholded
%   BF = regression_stats_results{c}.BF;    % statistic_image, type 'BF', 2*ln(BF10)
%
%   load(fullfile(resultsdir, 'svm_stats_results_<...>.mat'), 'svm_stats_results');
%   w  = svm_stats_results{c}.weight_obj;   % fmri_data: SVM weights, e.g. for 'topography'
%
%   load(fullfile(resultsdir, 'contrast_data_objects.mat'), 'DATA_OBJ_CON');
%   subj = DATA_OBJ_CON{c};                  % fmri_data_st, one image per participant:
%                                            % 'small_multiples', or 'image', k
%
%   load(fullfile(maskdir, '<model>_mask_<name>.mat'), 'mask', 'combined_atlas');
%                                            % from atlas_mask_tools: overlays and atlases
%
% LaBGAS maps are fmriprep 20 output: leave 'sourcespace' at its default.
