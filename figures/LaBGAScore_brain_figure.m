function [fig, o, info] = LaBGAScore_brain_figure(obj, varargin)
% LaBGAScore_brain_figure  Publication-quality brain figure from CANlab objects.
%
% Draws slice montages, cortical surfaces, 3-D cutaways and subcortical
% scenes, per-region close-ups, zoom insets, small multiples and cluster
% "topography" surfaces from the key CANlab object classes, alone or
% assembled into one multi-panel figure, and exports it at a fixed physical
% size and resolution that does not depend on the screen it runs on.
%
% It is the publication counterpart of plugin_set_figure_size, which sizes
% figures for publish() html reports by fitting them to the SCREEN. A paper
% figure needs the opposite: a size in mm, fonts in pt and a fixed dpi,
% identical whoever runs it, headless or not. Run it headless by default.
%
% USAGE:
% [fig, o, info] = LaBGAScore_brain_figure(obj, 'layout', LAYOUT, Name, Value, ...)
%
% LaBGAScore_brain_figure(t, 'layout', 'compact2', 'save', 'fig2a')
% LaBGAScore_brain_figure(t, 'layout', 'surfaces', 'threshold', {.05, 'fdr'})
% LaBGAScore_brain_figure(t, 'layout', 'slices', 'style', 'unthresholded_outline', 'sig', {.05, 'fdr'})
% LaBGAScore_brain_figure(BF, 'layout', 'slices')                    % Bayes-factor map
% LaBGAScore_brain_figure(w, 'layout', 'topography', 'roi', {'VStriatum'})
% LaBGAScore_brain_figure(t, 'layout', 'composite', 'panels', ...
%     {{'slices'}, {'surfaces'}, {'subcortex'}}, 'save', 'fig3', 'formats', {'png', 'pdf'})
%
% See LaBGAScore_brain_figure_example.m for one worked section per layout,
% and figures/README.md for the design decisions behind the defaults.
%
% -------------------------------------------------------------------------
% INPUT
% -------------------------------------------------------------------------
% obj   any of
%       - statistic_image      thresholded or not; its .sig is respected.
%                              .type 'BF' selects the Bayes-factor display.
%       - fmri_data / fmri_data_st / image_vector
%       - region               e.g. the region_objs_* saved by prep_3a
%       - atlas                drawn in unique colours, with a legend
%       - fmri_mask_image      drawn as one solid colour
%       - a .nii / .nii.gz path, read with fmri_data
%       - a cell array of the above, for 'multirow', 'small_multiples' and
%         'winner_take_all'
%
%       An object holding more than one image ERRORS unless 'image' says
%       which one to show (an index, or 'mean'). CANlab's own surface()
%       silently plots the mean and addblobs silently plots the first image;
%       neither is acceptable for a figure that goes into a paper.
%
% -------------------------------------------------------------------------
% LAYOUTS ('layout')
% -------------------------------------------------------------------------
% 'compact' 'compact2' 'compact3' 'full' 'full2' 'full hcp' 'allslices'
%                   CanlabCore's prepackaged montages (canlab_results_fmridisplay),
%                   with this function's colours, sizing and export.
% 'multirow'        one compact2 row per map; obj must be a cell array.
% 'slices'          your own slice rows, see 'slices'.
% 'regioncenters'   one close-up axis per region, see 'max_regions'.
% 'surfaces'        cortical surfaces, accurately projected, see 'surface'.
% 'cutaway'         3-D cutaway or slab, see 'cutaway'.
% 'subcortex'       3-D subcortical scene, see 'structures'. An atlas input
%                   is drawn as isosurfaces in unique colours instead.
% 'zoom'            slices with boxed, magnified insets, see 'zoom'.
% 'topography'      cluster surface plot, see 'topo_mode'.
% 'small_multiples' one identical panel per image, see 'panel'.
% 'composite'       several of the layouts above on one figure, see 'panels'.
%                   Montage and surface panels share ONE fmridisplay object,
%                   so they share one colour mapping and one colorbar.
%
% -------------------------------------------------------------------------
% OPTIONAL NAME-VALUE ARGUMENTS
% -------------------------------------------------------------------------
% COMPOSITE
% 'panels'      cell array of panel specs. Each spec is a cell
%               {layout, Name, Value, ...}; the name-value pairs override
%               the global options for that panel only, e.g.
%               {'slices', 'slices', {{'axial', -20:10:50}}} or
%               {'topography', 'map', svm_weights, 'roi', {'VStriatum'}}.
%               'map' gives a panel its own input object.
% 'grid'        [nrows ncols] for the panels, filled row by row
%               (default: one panel per row)
% 'heights'     relative row heights (default: each row's natural height)
% 'widths'      relative column widths (default: equal)
%
% INPUT HANDLING
% 'image'       which image of a multi-image object: an index, or 'mean'
% 'threshold'   cell passed to statistic_image/threshold before display,
%               e.g. {.05, 'fdr'}, {.001, 'unc', 'k', 10}. Default {}: the
%               object is shown as it is.
% 'map'         panel-level replacement for obj (composite panels only)
% 'underlay'    anatomical image under the slices: a canlab_get_underlay_image
%               keyword ('fmriprep20', 'icbm2009c_0.5mm', 'mni152_1mm',
%               'colin', 'keuken', ...) or a file. Default 'fmriprep20' =
%               fmriprep20_template.nii.gz, MNI152NLin2009cAsym - the space
%               fMRIPrep normalises LaBGAS data to, and CanlabCore's default.
% 'clip_mask'   mask (file or object) the map is restricted to before
%               display, e.g. which('brain_mask_fmriprep20_template_1000.nii')
%               or a gray-matter mask. Default []: shown as estimated. It
%               removes whole data voxels; a coarse voxel straddling the
%               edge is kept and still overhangs the 1 mm underlay.
%
% COLOUR
% 'colormode'   'auto' (default) | 'split' | 'ramp' | 'lut' | 'solid' |
%               'indexed' | 'categorical' | 'discrete' | 'bayes'.
%               'auto' picks: atlas -> indexed; fmri_mask_image -> solid;
%               statistic_image of type 'BF' -> bayes; 'categories' given ->
%               categorical; all values >= 0 -> ramp; otherwise split.
% 'cmaprange'   colour limits. 2 or 4 elements for split ([-max max] or
%               [-max -min min max]), 2 for ramp/lut. Default: from the
%               displayed values (99th percentile of |value| as the
%               maximum, the threshold as the minimum), recorded in
%               info.cmaprange. FIX THIS across figures that are compared.
% 'splitcolor'  {minneg maxneg minpos maxpos}, default CANlab hot/cool
%               {[0 0 1] [0 1 1] [1 .5 0] [1 1 0]}
% 'mincolor' / 'maxcolor'   ramp end colours (default red -> yellow)
% 'colormap'    for 'lut': a colormap name ('canlabCmap', 'parula', ...)
%               or an n x 3 matrix. Default canlabCmap.
% 'color'       for 'solid' (default [1 0 0])
% 'categories'  names for an integer-coded map (value k = category k), for
%               winner-take-all, zone or conjunction maps
% 'category_colors'  n x 3, one row per category / atlas region
% 'winner_take_all'  true with a cell array of maps: each voxel gets the
%               index of the map with the largest value, if that value
%               exceeds 'wta_threshold' (default 0)
% 'nlevels'     for 'discrete' (graded integer maps, e.g. n subjects):
%               number of levels (default: max value)
% 'atlas_level' colour an atlas by a coarser label level, e.g. 'labels_3':
%               parcels sharing a name there share a colour and one legend
%               entry (36 fine canlab2024 parcels of VStriatum, Amygdala and
%               insula become 10). Default '': one colour per parcel.
%
% BAYES FACTORS (colormode 'bayes')
% 'bf_scale'    scale the map is stored on: '2lnbf' (default; what
%               estimateBayesFactor and the LaBGAS templates store),
%               'lnbf' or 'bf10'. BF10 = exp(value/2) on the default scale,
%               NEVER exp(value) - see LaBGAScore/CLAUDE.md.
% 'bf_threshold' evidence threshold as BF10 (default 3, matching prep_3a's
%               |BF| > 3). Voxels with BF10 >= 3 are red (evidence for H1),
%               BF10 <= 1/3 blue (evidence for H0).
% 'bf_cap'      BF10 at which the colours saturate (default 100)
% 'bf_show'     'evidence' (default; only voxels beyond either threshold) or
%               'all' (weak evidence drawn in neutral grey)
%
% STYLE
% 'style'       'blobs' (default) or 'unthresholded_outline': the whole
%               map drawn pale/transparent, the significant voxels drawn at
%               full colour and outlined, as in Acil et al. 2026 Fig 2.
% 'sig'         significance for 'unthresholded_outline': a threshold cell
%               ({.05, 'fdr'}), or a thresholded statistic_image, region or
%               mask whose voxels count as significant (e.g. bootstrapped
%               SVM weights, which carry no p-values of their own)
% 'unthresholded_alpha'  opacity of the pale layer (default 0.35)
% 'outline_color' / 'outline_width'  (default black, 0.5 pt)
% 'thresholds'  multi-threshold layers, lenient first, e.g.
%               {{.01, 'unc'}, {.05, 'fdr'}} - each level gets its own
%               colours and a legend entry (Spisak et al. 2026 Fig 3)
% 'threshold_colors'  {pos1 neg1; pos2 neg2; ...}, one row per level
% 'overlay'     extra layers: {{obj2, 'outline', [0 0 0]}, {mask, 'solid',
%               [0 1 0]}}, drawn on top of the map
%
% SLICES
% 'slices'      rows of slices: a cell of rows, each row {orientation,
%               coords_mm} or a cell of such segments to mix orientations
%               in one row, e.g. {{{'sagittal', 0}, {'axial', -30:10:60}}}.
%               Default {{'sagittal', 0}, {'axial', -30:10:60}} in one row.
% 'slice_coords'  label each slice 'z = 10' (default true)
% 'lr_labels'   'L' / 'R' on axial and coronal rows (default true).
%               CANlab montages are in NEUROLOGICAL convention: left is left.
% 'slice_locator'  small midline sagittal with a line at each axial /
%               coronal slice (default false)
% 'brain_outline'  outline of the brain mask on every slice (default false)
%
% SURFACES
% 'surface'     'inflated' (default; HCP fsLR_32k) | 'flat' (fsLR_32k) |
%               'midthickness' (MNI152NLin2009cAsym mesh, no projection
%               needed because its vertices ARE in MNI space)
% 'views'       subset of {'left lateral', 'left medial', 'right medial',
%               'right lateral'} (default all four, in that order)
% 'projection'  true (default): MNIsurf volume-to-surface projection
%               (Wu et al. 2018) on fsLR meshes. Without it the default
%               rendering samples the volume at vertex positions that do not
%               lie in the MNI cortical ribbon and misses most of the signal.
% 'sourcespace' 'MNI152NLin2009cAsym' (default; fmriprep 20 output) or
%               'MNI152NLin6Asym' (SPM / FSL normalisation)
%
% 3-D
% 'cutaway'     addbrain cutaway keyword (default 'left_cutaway'); also
%               'right_cutaway', 'left_insula_slab', 'right_insula_slab',
%               'accumbens_slab', 'coronal_slabs_4', ...
% 'structures'  addbrain keywords composed into one subcortical scene
%               (default {'thalamus', 'caudate', 'putamen', 'amygdala', 'brainstem'})
% 'view3d'      [azimuth elevation] for cutaway / subcortex (default per keyword)
%
% REGION CLOSE-UPS
% 'max_regions' at most this many regions (default 20, the LaBGAS report
%               limit). More are reduced to the strongest by peak |value|,
%               with a warning; Inf shows all.
% 'orientation' slice orientation for region close-ups (default 'axial')
%
% LABELS
% 'labels'      false (default) | 'auto' (regions of the displayed map,
%               named from 'label_atlas') | n x 2 cell {name, [x y z]}.
%               Drawn beside each panel with dotted leader lines, in each
%               panel where the point is visible.
% 'label_atlas' atlas keyword (default 'canlab2024_fine_fmriprep20_2mm')
% 'label_field' label level the names come from (default 'labels_3'). In the
%               fine canlab2024 atlas: 'labels' 570 parcel names
%               (Ctx_p32pr_L), 'labels_2' 518, 'labels_3' 246
%               (midcingulate_L), 'labels_4' 83 regions (cingulate_ACC_mPFC_L);
%               'labels_5' is the source atlas of each parcel, not anatomy.
% 'max_labels'  at most this many automatic labels (default 8)
%
% ZOOM INSETS
% 'zoom'        what the insets zoom out of, one of:
%               - boxes: a cell of main slices, each {orientation, coord_mm,
%                 boxes, names}: boxes is k x 4 [u1 u2 v1 v2] in mm in the
%                 slice plane (axial x,y; coronal x,z; sagittal y,z), names
%                 a 1 x k cellstr. Kragel et al. 2023 Fig 1a.
%               - 'clusters': the 'max_zoom' strongest clusters of the
%                 (thresholded) map, each on a slice through its peak
%                 ('orientation'), boxed 5 mm beyond its extent and named
%                 from 'label_atlas'.
%               - ROI names in 'roi_atlas', e.g. {'VStriatum', 'Amygdala'}:
%                 a slice through each ROI's centre, the ROI outlined.
% 'zoom_inset'  'slice' (default): each inset is the boxed part of the slice,
%               magnified. 'topography': a 3-D surface of the UNTHRESHOLDED
%               values - the boxed plane for boxes and clusters, a spline
%               fit over the slice's two axes for ROIs - all insets on one
%               shared height / colour scale with their own colorbar.
% 'max_zoom'    at most this many clusters for 'zoom', 'clusters' (default 4)
%
% TOPOGRAPHY (cluster surface plots)
% 'topo_mode'   'roi' (default): ROI voxels projected onto two axes and
%               fitted with a thin-plate smoothing spline, with contours and
%               gradient arrows on the floor (Kragel et al. 2023 Fig 1c).
%               'slice': the raw values of one slice, cropped to a box
%               (generalises figures/cluster_surface_plots.m). Both read the
%               UNTHRESHOLDED values, even of a thresholded statistic_image:
%               a topography shows the values around a result.
% 'roi'         atlas | region | fmri_mask_image | cellstr of labels in
%               'roi_atlas' (default 'canlab2024_fine_fmriprep20_2mm')
% 'projection_axes'  two of 'x' 'absx' 'y' 'z' (default {'absx', 'y'};
%               absx folds the hemispheres, medial to lateral)
% 'hemisphere'  'both' (default) | 'left' | 'right'
% 'smoothing'   tpaps smoothing parameter p in [0 1]; default [] lets
%               tpaps choose, and the value used is in info.topography
% 'mask_distance'  hide the fitted surface further than this many mm from
%               any voxel (default 3); Inf shows the whole rectangle
% 'plane' / 'box'  for 'slice' mode: {orientation, coord_mm} and
%               [u1 u2 v1 v2] in mm, as for 'zoom'
% 'floor'       contours + gradient arrows below the surface (default: on
%               for 'topography', off for the smaller 'zoom' insets)
% 'inset'       3-D view of the ROI coloured by value (default true)
% 'zlabel'      height axis label (default 'value')
% 'topo_view'   [azimuth elevation] (default [-37.5 30])
%
% SMALL MULTIPLES
% 'panel'       what each image shows: {'axial', z}, {'coronal', y},
%               {'sagittal', x} or {'surface', 'left lateral'}
% 'ncols'       columns (default: ceil(sqrt(n)))
% 'titles'      cellstr, one per image (default: 'Image k')
% 'markers'     n x 3 mm coordinates marked with a ring on every panel
% 'marker_color' default [1 0 0]
%
% FIGURE AND EXPORT
% 'width'       'single' (89 mm) | 'onehalf' (120 mm) | 'double' (183 mm,
%               default) | a width in mm
% 'height'      mm (default: from the layout's natural aspect)
% 'fontsize'    pt (default 7) / 'fontname' (default 'Helvetica'; Arial is
%               not installed on the LaBGAS server, Liberation Sans is its
%               metric-compatible substitute)
% 'title'       figure title
% 'panel_letters'  'A' (default) | 'a' | false
% 'colorbar'    'auto' (default) | 'south' | 'east' | 'none'
% 'colorbar_label'  default from the map: 't', 'BF_{10}', 'weight', ...
% 'save'        file name without extension; '' (default) saves nothing
% 'formats'     cellstr from {'png', 'pdf', 'svg', 'tif'}, default {'png'}
% 'dpi'         default 600
% 'vector'      true writes pdf/svg text and lines as vectors (slow, ~4x,
%               and MATLAB warns it may misrender lit surfaces); false
%               (default) embeds the whole figure as an image at 'dpi'
% 'visible'     show the figure (default false)
% 'verbose'     default true
%
% -------------------------------------------------------------------------
% OUTPUT
% -------------------------------------------------------------------------
% fig   figure handle
% o     the fmridisplay object (empty for layouts that do not use one)
% info  what was drawn, for the figure legend and methods: colour mode and
%       limits, thresholds, regions and labels, the projection used, Bayes
%       factor ranges in BF10 units, topography fit parameters, files written
%       and any warnings
%
% -------------------------------------------------------------------------
% NOTES
% -------------------------------------------------------------------------
% - Every panel that shows the same map uses ONE canlab_colormap mapping,
%   and the colorbar is drawn from that same object, so it cannot disagree
%   with the rendering.
% - Surfaces are lit 3-D patches. In pdf/svg output they are always raster,
%   whatever 'vector' says.
% - Headless MATLAB has no OpenGL and renders with painters. Both renderers
%   were compared on the LaBGAS server; see figures/README.md.
% - Blobs past the cortical edge are usually resolution, not misregistration:
%   a 2-3 mm voxel drawn on the 1 mm underlay overhangs it. Data in another
%   space (e.g. older SPM/MNI152NLin6Asym data on the 2009cAsym default)
%   overhang more; pass a matching 'underlay', or 'clip_mask'.
%
% See also: canlab_results_fmridisplay, fmridisplay, addblobs, render_on_surface,
% addbrain, canlab_colormap, tpaps, plugin_set_figure_size
%
% -------------------------------------------------------------------------
%
% author: Lukas Van Oudenhove & Claude Opus 5.5
%
% date:   KU Leuven, November, 2026
%
% -------------------------------------------------------------------------
%
% LaBGAScore_brain_figure.m         v1.0
%
% last modified: 2026/11/26
%
%

%% PARSE OPTIONS
% -------------------------------------------------------------------------

p = build_parser();
parse(p, varargin{:});
P = p.Results;

if ~iscell(P.formats), P.formats = cellstr(P.formats); end
P.layout = lower(char(P.layout));
P.underlay = resolve_underlay(P.underlay);

info = struct('layout', P.layout, 'underlay', P.underlay, 'panels', {{}}, 'warnings', {{}}, 'files', {{}});
setappdata(groot, 'LaBGAScore_brain_figure_warnings', {});      % collected by lbf_warn

% fail early on the toolbox the topography layout needs, rather than after
% minutes of rendering
if any(strcmp(P.layout, {'topography'})) && isempty(which('tpaps')) && strcmp(P.topo_mode, 'roi')
    error('LaBGAScore_brain_figure:noCurveFitting', ...
        ['The ''roi'' topography fits a thin-plate smoothing spline with tpaps, ' ...
        'from the Curve Fitting Toolbox, which is not on the path. Use ''topo_mode'', ''slice'' instead.']);
end


%% PREPARE THE INPUT
% -------------------------------------------------------------------------

obj = prepare_input(obj, P);


%% PRESET MONTAGES
% -------------------------------------------------------------------------
% CanlabCore's own layouts open their own figure and place their own axes.
% They are kept as they are and only recoloured, resized and exported here;
% for full control over a layout use 'slices' / 'composite' instead.

presets = {'compact', 'compact2', 'compact3', 'full', 'full2', 'full hcp', 'full hcp inflated', 'allslices', 'multirow'};

if any(strcmp(P.layout, presets))
    [fig, o, info] = draw_preset(obj, P, info);
    info = finish_figure(fig, P, info);
    return
end


%% BUILD THE PANEL LIST
% -------------------------------------------------------------------------
% Every other layout is a composite with one panel, so all of them share the
% same layout, painting, labelling and export code.

if strcmp(P.layout, 'composite')
    if isempty(P.panels)
        error('LaBGAScore_brain_figure:noPanels', 'layout ''composite'' needs ''panels''.');
    end
    specs = P.panels;
elseif strcmp(P.layout, 'small_multiples')
    specs = small_multiple_specs(obj, P);
else
    specs = {{P.layout}};
end

panels = cell(1, numel(specs));
for i = 1:numel(specs)
    panels{i} = make_panel(specs{i}, obj, P, varargin);
end

% panel letters label panels of a multi-panel figure; a single panel gets none
% unless asked for
if numel(panels) == 1 && ~any(strcmpi(varargin(1:2:end), 'panel_letters'))
    panels{1}.Q.panel_letters = false;
end

% colour specification per distinct map. Panels showing the global map share
% one; small multiples share one computed over ALL their images.
[panels, info] = resolve_panel_colours(panels, obj, P, info);


%% FIGURE AND GRID
% -------------------------------------------------------------------------

[fig, G] = open_figure(panels, P);

% 'nocontroller': fmridisplay otherwise opens its interactive controller as a
% hidden uifigure for every object - in -nodisplay sessions too, where only
% -batch suppresses it - and they accumulate across calls
o = fmridisplay('overlay', P.underlay, 'nocontroller');
o.additional_info = 'LaBGAScore_brain_figure';
overlay = make_overlay_axes(fig);


%% DRAW EACH PANEL
% -------------------------------------------------------------------------

for i = 1:numel(panels)
    [panels{i}, o] = draw_panel(panels{i}, o, fig, G.rects{i}, P);
end


%% PAINT THE MAPS
% -------------------------------------------------------------------------

for i = 1:numel(panels)
    o = paint_panel(o, panels{i}, P);
end

% zoom limits are set AFTER painting, which may reset them
for i = 1:numel(panels)
    apply_axis_limits(panels{i});
end


%% ANNOTATE
% -------------------------------------------------------------------------

drawnow;                                   % camera and limits must be final before projecting points
uistack(overlay, 'top');

for i = 1:numel(panels)
    annotate_panel(panels{i}, overlay, fig, P);
end

info = draw_colour_keys(panels, overlay, G, fig, P, info);

if ~isempty(P.title)
    text(overlay, 0.5, 1 - G.mm2y(1.5), P.title, 'Units', 'normalized', ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
        'FontWeight', 'bold', 'FontSize', P.fontsize + 1);
end

info.panels = cellfun(@panel_summary, panels, 'UniformOutput', false);
info = finish_figure(fig, P, info);

end % main function



%% =======================================================================
%  OPTIONS
%  =======================================================================

function p = build_parser()

isnum = @(x) isnumeric(x) || islogical(x);
ischarlike = @(x) ischar(x) || isstring(x);

p = inputParser;
p.FunctionName = 'LaBGAScore_brain_figure';

% layout and composite
addParameter(p, 'layout', 'compact2', ischarlike);
addParameter(p, 'panels', {}, @iscell);
addParameter(p, 'grid', [], isnum);
addParameter(p, 'heights', [], isnum);
addParameter(p, 'widths', [], isnum);

% input handling
addParameter(p, 'image', [], @(x) isnumeric(x) || ischarlike(x));
addParameter(p, 'threshold', {}, @iscell);
addParameter(p, 'underlay', 'fmriprep20', ischarlike);
addParameter(p, 'clip_mask', [], @(x) isempty(x) || ischarlike(x) || isa(x, 'image_vector'));
addParameter(p, 'map', [], @(x) true);

% colour
addParameter(p, 'colormode', 'auto', ischarlike);
addParameter(p, 'cmaprange', [], isnum);
addParameter(p, 'splitcolor', {[0 0 1] [0 1 1] [1 .5 0] [1 1 0]}, @iscell);
addParameter(p, 'mincolor', [1 0 0], isnum);
addParameter(p, 'maxcolor', [1 1 0], isnum);
addParameter(p, 'colormap', [], @(x) isempty(x) || ischarlike(x) || isnumeric(x));
addParameter(p, 'color', [1 0 0], isnum);
addParameter(p, 'categories', {}, @iscell);
addParameter(p, 'category_colors', [], isnum);
addParameter(p, 'winner_take_all', false, isnum);
addParameter(p, 'wta_threshold', 0, isnum);
addParameter(p, 'nlevels', [], isnum);
addParameter(p, 'atlas_level', '', ischarlike);

% Bayes factors
addParameter(p, 'bf_scale', '2lnbf', ischarlike);
addParameter(p, 'bf_threshold', 3, isnum);
addParameter(p, 'bf_cap', 100, isnum);
addParameter(p, 'bf_show', 'evidence', ischarlike);

% style
addParameter(p, 'style', 'blobs', ischarlike);
addParameter(p, 'sig', {}, @(x) true);
addParameter(p, 'unthresholded_alpha', 0.35, isnum);
addParameter(p, 'outline_color', [0 0 0], isnum);
addParameter(p, 'outline_width', 0.5, isnum);
addParameter(p, 'thresholds', {}, @iscell);
addParameter(p, 'threshold_colors', {}, @iscell);
addParameter(p, 'overlay', {}, @iscell);

% slices
addParameter(p, 'slices', {}, @iscell);
addParameter(p, 'slice_coords', true, isnum);
addParameter(p, 'lr_labels', true, isnum);
addParameter(p, 'slice_locator', false, isnum);
addParameter(p, 'brain_outline', false, isnum);

% surfaces
addParameter(p, 'surface', 'inflated', ischarlike);
addParameter(p, 'views', {'left lateral', 'left medial', 'right medial', 'right lateral'}, @iscell);
addParameter(p, 'projection', true, isnum);
addParameter(p, 'sourcespace', 'MNI152NLin2009cAsym', ischarlike);

% 3-D
addParameter(p, 'cutaway', 'left_cutaway', ischarlike);
addParameter(p, 'structures', {'thalamus', 'caudate', 'putamen', 'amygdala', 'brainstem'}, @iscell);
addParameter(p, 'view3d', [], isnum);

% region close-ups
addParameter(p, 'max_regions', 20, isnum);
addParameter(p, 'orientation', 'axial', ischarlike);

% labels
addParameter(p, 'labels', false, @(x) islogical(x) || ischarlike(x) || iscell(x));
addParameter(p, 'label_atlas', 'canlab2024_fine_fmriprep20_2mm', ischarlike);
addParameter(p, 'label_field', 'labels_3', ischarlike);
addParameter(p, 'max_labels', 8, isnum);

% zoom
addParameter(p, 'zoom', {}, @(x) iscell(x) || ischarlike(x));
addParameter(p, 'zoom_inset', 'slice', ischarlike);
addParameter(p, 'max_zoom', 4, isnum);

% topography
addParameter(p, 'topo_mode', 'roi', ischarlike);
addParameter(p, 'roi', [], @(x) true);
addParameter(p, 'roi_atlas', 'canlab2024_fine_fmriprep20_2mm', ischarlike);
addParameter(p, 'projection_axes', {'absx', 'y'}, @iscell);
addParameter(p, 'hemisphere', 'both', ischarlike);
addParameter(p, 'smoothing', [], isnum);
addParameter(p, 'mask_distance', 3, isnum);
addParameter(p, 'plane', {}, @iscell);
addParameter(p, 'box', [], isnum);
addParameter(p, 'floor', [], isnum);
addParameter(p, 'inset', true, isnum);
addParameter(p, 'zlabel', 'value', ischarlike);
addParameter(p, 'topo_view', [-37.5 30], isnum);

% small multiples
addParameter(p, 'panel', {'axial', 0}, @iscell);
addParameter(p, 'ncols', [], isnum);
addParameter(p, 'titles', {}, @iscell);
addParameter(p, 'markers', [], isnum);
addParameter(p, 'marker_color', [1 0 0], isnum);
addParameter(p, 'image_index', [], isnum);          % internal: small multiple panel k

% figure and export
addParameter(p, 'width', 'double', @(x) isnumeric(x) || ischarlike(x));
addParameter(p, 'height', [], isnum);
addParameter(p, 'fontsize', 7, isnum);
addParameter(p, 'fontname', 'Helvetica', ischarlike);
addParameter(p, 'title', '', ischarlike);
addParameter(p, 'panel_title', '', ischarlike);
addParameter(p, 'panel_letters', 'A', @(x) ischarlike(x) || islogical(x));
addParameter(p, 'colorbar', 'auto', ischarlike);
addParameter(p, 'colorbar_label', '', ischarlike);
addParameter(p, 'save', '', ischarlike);
addParameter(p, 'formats', {'png'}, @(x) iscell(x) || ischarlike(x));
addParameter(p, 'dpi', 600, isnum);
addParameter(p, 'vector', false, isnum);
addParameter(p, 'visible', false, isnum);
addParameter(p, 'verbose', true, isnum);

end



%% =======================================================================
%  INPUT
%  =======================================================================

function obj = prepare_input(obj, P)
% Read paths, resolve cell inputs and pick the image to show.

if ischar(obj) || isstring(obj)
    obj = read_image(obj, P);

elseif iscell(obj)
    for i = 1:numel(obj)
        if ischar(obj{i}) || isstring(obj{i}), obj{i} = read_image(obj{i}, P); end
        obj{i} = clip_to_mask(obj{i}, P.clip_mask);
    end
    if P.winner_take_all
        obj = winner_take_all(obj, P);
    elseif ~any(strcmp(P.layout, {'multirow', 'small_multiples', 'composite'}))
        error('LaBGAScore_brain_figure:cellInput', ...
            'A cell array of maps needs layout ''multirow'' or ''small_multiples'', or ''winner_take_all'', true.');
    end
    return
end

obj = select_image(obj, P);
if isa(obj, 'atlas') && ~isempty(P.atlas_level)
    obj = merge_atlas_level(obj, P.atlas_level);
end
obj = clip_to_mask(obj, P.clip_mask);

end


function f = resolve_underlay(u)
% Keyword or file -> full path, so a missing underlay fails here and not
% halfway through drawing.

u = char(u);
if exist(u, 'file') == 2
    f = u;
else
    f = canlab_get_underlay_image(u);
end
if isempty(f)
    error('LaBGAScore_brain_figure:noUnderlay', 'Underlay not found: %s', u);
end

end


function obj = clip_to_mask(obj, mask)
% apply_mask also drops .sig/.p of a statistic_image with the voxels.

if isempty(mask), return, end
if isa(obj, 'region')
    lbf_warn('LaBGAScore_brain_figure:clipRegion', '''clip_mask'' is ignored for region objects; pass the image instead.');
    return
end
if ischar(mask) || isstring(mask)
    mask = fmri_mask_image(char(mask), 'noverbose');
end
obj = apply_mask(obj, mask);

end


function out = merge_atlas_level(atl, level)
% One parcel per distinct name at a coarser label level.

if ~isprop(atl, level) || isempty(atl.(level))
    error('LaBGAScore_brain_figure:atlasLevel', 'Atlas has no label level ''%s''.', level);
end
[names, ~, g] = unique(atl.(level), 'stable');
idx = round(double(atl.dat));
merged = zeros(size(idx));
nz = idx > 0 & idx <= numel(g);
merged(nz) = g(idx(nz));
out = atl;
out.dat = merged;
out.labels = names(:)';
out.probability_maps = [];
for f = {'labels_2', 'labels_3', 'labels_4', 'labels_5', 'label_descriptions'}
    if isprop(out, f{1}), out.(f{1}) = {}; end
end

end


function obj = read_image(fname, P)

fname = char(fname);
if ~exist(fname, 'file')
    error('LaBGAScore_brain_figure:noFile', 'Image file not found: %s', fname);
end
obj = fmri_data(fname, 'noverbose');
obj = select_image(obj, P);

end


function obj = select_image(obj, P)
% A figure shows ONE image. Refuse to guess which.

if ~isa(obj, 'image_vector') || isa(obj, 'atlas') || strcmp(P.layout, 'small_multiples')
    return
end

nimg = size(obj.dat, 2);
if nimg <= 1, return, end

if isempty(P.image)
    error('LaBGAScore_brain_figure:multiImage', ...
        ['The input holds %d images. Say which one to show with ''image'', k (or ''image'', ''mean''), ' ...
        'or use layout ''small_multiples'' to show them all.'], nimg);
end

if ischar(P.image) || isstring(P.image)
    if ~strcmpi(P.image, 'mean')
        error('LaBGAScore_brain_figure:badImage', '''image'' must be an index or ''mean''.');
    end
    obj = mean(obj);
else
    obj = pick_image(obj, P.image);
end

end


function out = pick_image(obj, k)
% The images k of a multi-image object. Not get_wh_image: on an fmri_data_st
% (DATA_OBJ_CON in the LaBGAS templates) it returned all 21 images when asked
% for images 1:9.

out = obj;
n = size(obj.dat, 2);
out.dat = obj.dat(:, k);
if isa(obj, 'statistic_image')
    for f = {'p', 'ste', 'sig'}
        if size(obj.(f{1}), 2) == n, out.(f{1}) = obj.(f{1})(:, k); end
    end
end
names = obj.image_names;
if iscell(names) && numel(names) == n, out.image_names = names(k);
elseif ischar(names) && size(names, 1) == n, out.image_names = names(k, :); end
if ~isempty(obj.fullpath) && size(obj.fullpath, 1) == n, out.fullpath = obj.fullpath(k, :); end
if numel(obj.removed_images) == n, out.removed_images = obj.removed_images(k); end

end


function wta = winner_take_all(maps, P)
% Each voxel gets the index of the map with the largest value, where that
% value exceeds the threshold; 0 elsewhere.

base = maps{1};
for i = 2:numel(maps)
    maps{i} = resample_space(maps{i}, base);
end
dat = zeros(size(base.dat, 1), numel(maps));
for i = 1:numel(maps)
    d = double(maps{i}.dat(:, 1));
    if isa(maps{i}, 'statistic_image') && ~isempty(maps{i}.sig)
        d(~maps{i}.sig(:, 1)) = 0;
    end
    dat(:, i) = d;
end
[mx, idx] = max(dat, [], 2);
idx(mx <= P.wta_threshold) = 0;

wta = fmri_data(base, 'noverbose');
wta.dat = double(idx);
wta.image_names = 'winner_take_all';

end



%% =======================================================================
%  PANELS
%  =======================================================================

function specs = small_multiple_specs(obj, P)
% One panel per image; they share a colour specification computed over all.
% 'image', [list] picks which images (default all).

if iscell(obj)
    idx = 1:numel(obj);
else
    idx = 1:size(obj.dat, 2);
end
if isnumeric(P.image) && ~isempty(P.image), idx = P.image(:)'; end
specs = cell(1, numel(idx));
for j = 1:numel(idx)
    k = idx(j);
    title_k = sprintf('Image %d', k);
    if numel(P.titles) >= j, title_k = P.titles{j}; end
    % slice coordinates would repeat on every panel; L / R only on the first
    specs{j} = {'small_multiple', 'image_index', k, 'panel_title', title_k, 'panel_letters', false, ...
        'slice_coords', false, 'lr_labels', j == 1};
end

end


function pn = make_panel(spec, obj, ~, global_args)
% A panel = its layout type plus the options in force for it. Panel-level
% name-value pairs are appended to the global ones, and inputParser keeps the
% LAST value of a repeated option, so they override.

if ischar(spec), spec = {spec}; end
type = lower(spec{1});

pp = build_parser();
parse(pp, global_args{:}, spec{2:end});
Q = pp.Results;
Q.layout = type;

if ~isempty(Q.map)
    map = prepare_input(Q.map, Q);
    own_map = true;
else
    map = obj;
    own_map = false;
end

if strcmp(type, 'small_multiple')
    if iscell(obj)
        map = select_image(obj{Q.image_index}, setfield(Q, 'image', 1)); %#ok<SFLD>
    else
        map = pick_image(obj, Q.image_index);
    end
end

pn = struct('type', type, 'Q', Q, 'map', map, 'own_map', own_map, ...
    'C', [], 'layers', {{}}, 'montages', [], 'surf_fslr', [], 'surf_mni', [], ...
    'axes', gobjects(0), 'slice_info', struct('ax', {}, 'orientation', {}, 'coord', {}), ...
    'limits', struct('ax', {}, 'xlim', {}, 'ylim', {}), 'view_axes', gobjects(0), ...
    'regions', [], 'labels', {{}}, 'topography', [], 'rect', [], 'key_rect', [], 'label_band', [], ...
    'titles', struct('ax', {}, 'lines', {}));

end


function s = panel_summary(pn)

s = struct('type', pn.type, 'colormode', '', 'cmaprange', []);
if ~isempty(pn.C)
    s.colormode = pn.C.mode;
    s.cmaprange = pn.C.range;
end
if ~isempty(pn.topography), s.topography = pn.topography; end
if ~isempty(pn.regions), s.n_regions = numel(pn.regions); end
if ~isempty(pn.labels), s.labels = pn.labels; end

end



%% =======================================================================
%  COLOUR
%  =======================================================================

function [panels, info] = resolve_panel_colours(panels, obj, P, info)
% One colour specification per distinct map, so panels showing the same map
% agree exactly. Small multiples share one over ALL images.

global_C = [];
if strcmp(P.layout, 'small_multiples')
    pool = obj;
    if ~iscell(obj) && isnumeric(P.image) && ~isempty(P.image), pool = pick_image(obj, P.image); end
    global_C = resolve_colour(pool, P);
end

for i = 1:numel(panels)
    pn = panels{i};
    if any(strcmp(pn.type, {'topography'}))
        continue                                       % topography has its own colour handling
    end
    if strcmp(pn.type, 'small_multiple')
        pn.C = global_C;
    elseif pn.own_map || ~isempty(setdiff_colour_opts(pn.Q, P))
        pn.C = resolve_colour(pn.map, pn.Q);
    else
        if isempty(global_C), global_C = resolve_colour(obj, P); end
        pn.C = global_C;
    end
    pn.layers = build_layers(pn.map, pn.C, pn.Q);
    panels{i} = pn;
end

if ~isempty(global_C)
    info.colormode = global_C.mode;
    info.cmaprange = global_C.range;
    if isfield(global_C, 'bayes'), info.bayes = global_C.bayes; end
end

end


function d = setdiff_colour_opts(Q, P)
% Names of colour options a panel overrides, if any.

names = {'colormode', 'cmaprange', 'splitcolor', 'mincolor', 'maxcolor', 'colormap', 'color', ...
    'categories', 'category_colors', 'threshold', 'style', 'sig', 'thresholds', 'bf_threshold', 'bf_cap', 'bf_show'};
d = {};
for i = 1:numel(names)
    if ~isequal(Q.(names{i}), P.(names{i})), d{end+1} = names{i}; end %#ok<AGROW>
end

end


function C = resolve_colour(map, Q)
% Work out the colour mode, the addblobs arguments and a canlab_colormap
% object to draw the colorbar from.

if iscell(map)                                   % small multiples / multirow: pool the values
    vals = [];
    for i = 1:numel(map), vals = [vals; shown_values(apply_threshold(map{i}, Q))]; end %#ok<AGROW>
    probe = map{1};
else
    probe = apply_threshold(map, Q);
    vals = shown_values(probe);
    if size(probe.dat, 2) > 1 && ~isa(probe, 'atlas')
        vals = shown_values_all(probe);
    end
end

mode = lower(char(Q.colormode));
if strcmp(mode, 'auto')
    mode = auto_mode(probe, vals, Q);
end

C = struct('mode', mode, 'args', {{}}, 'cmap', [], 'range', [], 'ticks', [], 'ticklabels', {{}}, ...
    'label', char(Q.colorbar_label), 'legend', struct('name', {}, 'color', {}), 'thresholded', false);

if isempty(C.label), C.label = default_value_label(probe, mode); end

if isempty(vals), vals = [-1; 1]; end
if strcmpi(Q.style, 'unthresholded_outline') && ~iscell(map) && ~isempty(Q.sig)
    % the pale layer is every voxel, but the colours must span the SIGNIFICANT
    % ones, or all of them saturate at the end of the scale
    [~, sigmap] = split_significance(map, Q);
    sv = shown_values(sigmap);
    if ~isempty(sv), vals = [vals; sign(sv) .* max(abs(sv))]; end
end
absv = abs(vals(vals ~= 0));
if isempty(absv), absv = 1; end
vmax = robust_max(absv);
thresholded = is_thresholded(probe) || ~isempty(Q.threshold);
C.thresholded = thresholded;

switch mode

    case 'split'
        if ~isempty(Q.cmaprange)
            r = Q.cmaprange(:)';
            if numel(r) == 2, r = [min(r) 0 0 max(r)]; end
        else
            vmin = 0;
            if thresholded, vmin = min(absv); end
            r = [-vmax -vmin vmin vmax];
        end
        sc = Q.splitcolor;
        C.args = {'splitcolor', sc, 'cmaprange', r};
        C.cmap = canlab_colormap.split(sc{1}, sc{2}, sc{3}, sc{4}, r);
        C.range = r;
        if r(3) > 0
            C.ticks = [r(1) r(2) r(3) r(4)];
        else
            C.ticks = [r(1) 0 r(4)];
        end

    case 'ramp'
        r = Q.cmaprange;
        if isempty(r), r = [min(vals) robust_max(vals)]; end
        C.args = {'mincolor', Q.mincolor, 'maxcolor', Q.maxcolor, 'cmaprange', r};
        C.cmap = canlab_colormap.single(Q.mincolor, Q.maxcolor, r);
        C.range = r;
        C.ticks = r;

    case 'lut'
        cm = get_colormap(Q.colormap);
        r = Q.cmaprange;
        if isempty(r)
            if any(vals < 0), r = [-vmax vmax]; else, r = [min(vals) robust_max(vals)]; end
        end
        C.args = {'colormap', cm, 'cmaprange', r};
        C.cmap = canlab_colormap.continuous(cm, r);
        C.range = r;
        C.ticks = unique([r(1) 0 r(2)]);
        if r(1) >= 0, C.ticks = r; end

    case 'solid'
        C.args = {'color', Q.color};
        C.cmap = canlab_colormap.solid(Q.color);
        C.range = [];
        C.legend(1).name = C.label;
        C.legend(1).color = Q.color;

    case {'indexed', 'categorical'}
        if isa(probe, 'atlas')
            % colour by the index image, which is what gets drawn; .labels can
            % list parcels that .dat no longer holds (e.g. a merged atlas whose
            % parcels all carry index 1)
            names = probe.labels;
            n = num_regions(probe);
            present = unique(round(double(probe.dat(probe.dat ~= 0))));
            if numel(present) ~= n
                lbf_warn('LaBGAScore_brain_figure:atlasIndex', ...
                    ['Atlas lists %d labels but its index image (.dat) holds %d distinct parcels; ' ...
                    'drawing .dat, and listing only those parcels in the legend.'], n, numel(present));
            end
        else
            present = [];
            n = max(round(vals));
            names = Q.categories;
            if numel(names) < n
                names = [names, arrayfun(@(k) sprintf('%d', k), numel(names)+1:n, 'UniformOutput', false)];
            end
        end
        cmap = Q.category_colors;
        if isempty(cmap), cmap = default_category_colors(n, isa(probe, 'atlas')); end
        if size(cmap, 1) < n
            error('LaBGAScore_brain_figure:fewColors', '''category_colors'' has %d rows for %d categories.', size(cmap, 1), n);
        end
        cmap = cmap(1:n, :);
        C.args = {'indexmap', cmap, 'interp', 'nearest'};      % indexes must not be interpolated
        C.cmap = canlab_colormap.indexed(cmap);
        C.range = [1 n];
        shown = 1:n;
        if ~isempty(present), shown = intersect(shown, present(:)'); end
        if numel(shown) <= 24                            % beyond that a legend is unreadable
            for k = shown
                C.legend(end+1).name = strrep(names{k}, '_', ' ');
                C.legend(end).color = cmap(k, :);
            end
        end

    case 'discrete'
        n = Q.nlevels;
        if isempty(n), n = max(round(vals)); end
        cmap = sequential_steps(n);
        C.args = {'indexmap', cmap, 'interp', 'nearest'};
        C.cmap = canlab_colormap.indexed(cmap);
        C.range = [1 n];
        C.ticks = 1:n;

    case 'bayes'
        C = bayes_colour(C, vals, Q);

    otherwise
        error('LaBGAScore_brain_figure:badColormode', 'Unknown colormode ''%s''.', mode);
end

C.ticklabels = arrayfun(@num_label, C.ticks, 'UniformOutput', false);
if strcmp(mode, 'bayes'), C.ticklabels = C.bayes.ticklabels; end

end


function mode = auto_mode(probe, vals, Q)

if isa(probe, 'atlas')
    mode = 'indexed';
elseif isa(probe, 'fmri_mask_image')
    mode = 'solid';
elseif isa(probe, 'statistic_image') && strcmpi(probe.type, 'BF')
    mode = 'bayes';
elseif ~isempty(Q.categories) || Q.winner_take_all
    mode = 'categorical';
elseif ~isempty(vals) && all(vals >= 0)
    mode = 'ramp';
else
    mode = 'split';
end

end


function C = bayes_colour(C, vals, Q)
% Diverging map anchored at the evidence thresholds (Kwon et al. 2026 Fig 2a):
% blue = evidence for H0, grey = weak, red = evidence for H1, on the 2*ln(BF10)
% scale the LaBGAS templates store. The colorbar is labelled in BF10.

v = to_2lnbf(vals, Q.bf_scale);
T = 2 * log(Q.bf_threshold);
cap = 2 * log(Q.bf_cap);

lo = max(min(v), -cap);
hi = min(max(v), cap);
if hi <= lo, hi = lo + 1; end

n = 256;
x = linspace(lo, hi, n)';
lut = repmat([.85 .85 .85], n, 1);                       % weak evidence: neutral grey

dark_blue = [0.03 0.19 0.42]; light_blue = [0.62 0.79 0.88];
light_red = [0.99 0.73 0.63]; dark_red = [0.65 0.06 0.08];

key_lut = lut;
if strcmpi(Q.bf_show, 'evidence')
    key_lut(x > -T & x < T, :) = 1;                    % not drawn on the brain: white on the key
end

h0 = x <= -T;
if any(h0)
    w = (x(h0) - (-T)) ./ (-cap - (-T));                  % 0 at threshold, 1 at cap
    lut(h0, :) = (1 - min(w, 1)) .* light_blue + min(w, 1) .* dark_blue;
end
h1 = x >= T;
if any(h1)
    w = (x(h1) - T) ./ (cap - T);
    lut(h1, :) = (1 - min(w, 1)) .* light_red + min(w, 1) .* dark_red;
end

key_lut(h0 | h1, :) = lut(h0 | h1, :);
C.args = {'colormap', lut, 'cmaprange', [lo hi]};
C.cmap = canlab_colormap.continuous(key_lut, [lo hi]);
C.range = [lo hi];
C.label = 'BF_{10}';

bf_ticks = [1/100 1/30 1/10 1/3 1 3 10 30 100 300 1000];
xt = 2 * log(bf_ticks);
keep = xt >= lo - 1e-9 & xt <= hi + 1e-9;
C.ticks = xt(keep);
labs = {'1/100', '1/30', '1/10', '1/3', '1', '3', '10', '30', '100', '300', '1000'};
labs = labs(keep);
if max(v) > cap && ~isempty(labs), labs{end} = ['\geq' labs{end}]; end

C.bayes = struct('scale_stored', Q.bf_scale, 'threshold_bf10', Q.bf_threshold, 'cap_bf10', Q.bf_cap, ...
    'min_bf10', exp(min(v) / 2), 'max_bf10', exp(max(v) / 2), ...
    'n_voxels_h1', sum(v >= T), 'n_voxels_h0', sum(v <= -T), 'ticklabels', {labs}, 'note', '');

if C.bayes.min_bf10 > 1/10
    C.bayes.note = sprintf(['No voxel reaches BF10 < 1/10, strong evidence for H0 (minimum BF10 = %.3f). ' ...
        'With the JZS prior that is a floor set by the sample size - below roughly n = 100 it is unreachable ' ...
        'however null the data - not a weak result: see "Bayes factors" in LaBGAScore/CLAUDE.md.'], ...
        C.bayes.min_bf10);
end

end


function v = to_2lnbf(v, scale)

switch lower(char(scale))
    case '2lnbf', v = double(v);
    case 'lnbf',  v = 2 * double(v);
    case 'bf10',  v = 2 * log(double(v));
    otherwise, error('LaBGAScore_brain_figure:bfScale', '''bf_scale'' must be ''2lnbf'', ''lnbf'' or ''bf10''.');
end

end


function cm = get_colormap(spec)

if isempty(spec), spec = 'canlabCmap'; end
if ischar(spec) || isstring(spec)
    spec = char(spec);
    if strcmp(spec, 'canlabCmap')
        cm = canlabCmap;
    else
        cm = feval(spec, 256);
    end
else
    cm = spec;
end

end


function cmap = default_category_colors(n, isatlas)
% Distinct, colour-blind-aware categories; atlases fall back to CANlab's
% scn_standard_colors when there are more regions than a palette holds.

pal = [0.49 0.18 0.56; 0.85 0.22 0.20; 0.98 0.80 0.10; 0.16 0.29 0.70; ...
    0.30 0.75 0.40; 0.95 0.55 0.15; 0.25 0.75 0.85; 0.85 0.40 0.75; ...
    0.55 0.35 0.20; 0.45 0.45 0.45];
if n <= size(pal, 1) && ~isatlas
    cmap = pal(1:n, :);
else
    c = scn_standard_colors(n);
    cmap = cat(1, c{:});
end

end


function cmap = sequential_steps(n)
% Light cream to dark red, one step per level (Sun et al. 2026 Fig 6C).

lo = [0.99 0.96 0.86]; hi = [0.45 0.04 0.10];
w = linspace(0, 1, n)';
cmap = (1 - w) .* lo + w .* hi;

end


function s = default_value_label(probe, mode)

if strcmp(mode, 'bayes'), s = 'BF_{10}'; return, end
s = 'value';
if isa(probe, 'statistic_image') && ~isempty(probe.type)
    s = probe.type;
    if strcmpi(s, 'generic'), s = 'value'; end
    if strcmp(s, 'T'), s = '\itt'; end                 % statistic_image stores 'T'
elseif isa(probe, 'fmri_data') && ~isempty(probe.image_names)
    names = probe.image_names;
    if iscell(names), n = names{1}; else, n = names(1, :); end
    if contains(lower(char(n)), 'weight'), s = 'weight'; end
end

end


function s = num_label(x)

if x == 0, s = '0'; return, end
s = sprintf('%.3g', x);

end


function m = robust_max(a)
% 99th percentile, so a handful of extreme voxels does not wash out the scale.

a = a(isfinite(a));
if isempty(a), m = 1; return, end
m = prctile(a, 99);
if m <= 0, m = max(a); end
if m <= 0, m = 1; end

end



%% =======================================================================
%  LAYERS
%  =======================================================================

function layers = build_layers(map, C, Q)
% The ordered list of things addblobs paints for one panel.

layers = {};
if iscell(map) || any(strcmp(Q.layout, {'topography'}))
    return
end

if ~isempty(Q.thresholds)
    layers = multi_threshold_layers(map, Q);

elseif strcmpi(Q.style, 'unthresholded_outline')
    [unthr, sig] = split_significance(map, Q);
    % the pale layer leaves the significant voxels out, so the opaque layer
    % does not depend on being drawn on top of it
    in_sig = double(sig.dat(:, 1)) ~= 0;
    if isa(sig, 'statistic_image') && ~isempty(sig.sig), in_sig = in_sig & logical(sig.sig(:, 1)); end
    unthr.dat(in_sig, 1) = 0;
    if isa(unthr, 'statistic_image'), unthr.sig(in_sig, 1) = 0; end
    layers{end+1} = layer(unthr, [C.args, {'trans', 'transvalue', Q.unthresholded_alpha}], false);
    layers{end+1} = layer(sig, C.args, false);
    layers{end+1} = layer(sig, {'outline', 'color', Q.outline_color, 'linewidth', Q.outline_width}, true);

elseif strcmp(C.mode, 'bayes')
    layers{end+1} = layer(bayes_display_map(map, Q), C.args, false);

else
    layers{end+1} = layer(apply_threshold(map, Q), C.args, false);
end

for i = 1:numel(Q.overlay)
    ov = Q.overlay{i};
    src = ov{1}; kind = 'outline'; col = [0 0 0];
    if numel(ov) >= 2, kind = ov{2}; end
    if numel(ov) >= 3, col = ov{3}; end
    switch kind
        case 'outline', layers{end+1} = layer(src, {'outline', 'color', col, 'linewidth', Q.outline_width}, true); %#ok<AGROW>
        case 'solid',   layers{end+1} = layer(src, {'color', col}, false); %#ok<AGROW>
        case 'trans',   layers{end+1} = layer(src, {'color', col, 'trans'}, false); %#ok<AGROW>
        otherwise, error('LaBGAScore_brain_figure:badOverlay', 'overlay kind must be outline, solid or trans.');
    end
end

if Q.brain_outline
    layers{end+1} = layer(brain_mask(), {'outline', 'color', [0 0 0], 'linewidth', Q.outline_width}, true);
end

end


function L = layer(src, args, montage_only)

L = struct('src', src, 'args', {args}, 'montage_only', montage_only);

end


function layers = multi_threshold_layers(map, Q)
% Lenient level first, so stricter levels are painted on top of it.

if ~isa(map, 'statistic_image')
    error('LaBGAScore_brain_figure:multiThreshold', '''thresholds'' needs a statistic_image (it has p-values to threshold).');
end

cols = Q.threshold_colors;
if isempty(cols)
    cols = {[0.84 0.19 0.15], [0.17 0.33 0.85]; ...
        [1.00 0.85 0.00], [0.20 0.85 0.35]; ...
        [1.00 1.00 0.70], [0.70 1.00 1.00]};
end

% Each level shows only the voxels that do NOT also pass the next, stricter
% level. Stacking overlapping layers and relying on draw order does not work:
% CanlabCore never restacks blobs, and the stricter layer came out hidden
% under the lenient one (FDR 2190 voxels of p < .01's 3511 in emotionreg,
% almost none of them visible).
tk = cell(1, numel(Q.thresholds));
for k = 1:numel(Q.thresholds)
    spec = Q.thresholds{k};
    tk{k} = threshold(map, spec{:}, 'noverbose');
end
for k = 1:numel(tk) - 1
    stricter = logical(tk{k + 1}.sig(:, 1));
    tk{k}.sig(stricter, 1) = 0;
end

layers = {};
for k = 1:numel(Q.thresholds)
    spec = Q.thresholds{k};
    [pos, neg] = split_sign(tk{k});
    layers{end+1} = layer(pos, {'color', cols{k, 1}}, false); %#ok<AGROW>
    layers{end+1} = layer(neg, {'color', cols{k, 2}}, false); %#ok<AGROW>
    layers{end}.legend = {['Increase, ' threshold_text(spec)], cols{k, 1}; ['Decrease, ' threshold_text(spec)], cols{k, 2}};
end

end


function s = threshold_text(spec)

pval = spec{1};
kind = 'unc';
if numel(spec) >= 2, kind = lower(spec{2}); end
switch kind
    case 'fdr', s = sprintf('FDR q < %s', strip_zero(pval));
    case 'fwe', s = sprintf('FWE p < %s', strip_zero(pval));
    case {'bfr', 'bonferroni'}, s = sprintf('Bonferroni p < %s', strip_zero(pval));
    otherwise,  s = sprintf('p < %s unc.', strip_zero(pval));
end

end


function s = strip_zero(x)

s = sprintf('%g', x);
if strncmp(s, '0.', 2), s = s(2:end); end

end


function [pos, neg] = split_sign(t)

pos = t; neg = t;
d = double(t.dat(:, 1));
pos.dat(:, 1) = max(d, 0);
neg.dat(:, 1) = min(d, 0);
if isa(t, 'statistic_image') && ~isempty(t.sig)
    pos.sig(:, 1) = t.sig(:, 1) & d > 0;
    neg.sig(:, 1) = t.sig(:, 1) & d < 0;
end

end


function [unthr, sig] = split_significance(map, Q)
% The pale layer is every non-zero voxel; the significant layer is the subset
% that passes Q.sig.

unthr = map;
if isa(unthr, 'statistic_image'), unthr.sig = true(size(unthr.dat)); end

s = Q.sig;
if isempty(s)
    error('LaBGAScore_brain_figure:noSig', 'style ''unthresholded_outline'' needs ''sig''.');
end

if iscell(s)
    if ~isa(map, 'statistic_image')
        error('LaBGAScore_brain_figure:sigNeedsStat', ...
            ['A threshold cell for ''sig'' needs a statistic_image. For a map without p-values ' ...
            '(e.g. SVM weights) pass the significant voxels instead: a thresholded statistic_image, region or mask.']);
    end
    sig = threshold(unthr, s{:}, 'noverbose');
    return
end

% a mask-like object: keep the map's values where it is non-zero / significant
if isa(s, 'region'), s = region2imagevec(s); end
s = resample_space(s, map, 'nearest');
keep = double(s.dat(:, 1)) ~= 0;
if isa(s, 'statistic_image') && ~isempty(s.sig), keep = keep & s.sig(:, 1); end
sig = map;
sig.dat(~keep, 1) = 0;
if isa(sig, 'statistic_image'), sig.sig = keep; end

end


function out = bayes_display_map(map, Q)
% Convert to the 2*ln(BF10) scale; in 'evidence' mode keep only voxels beyond
% either threshold.

out = map;
v = to_2lnbf(map.dat(:, 1), Q.bf_scale);
if strcmpi(Q.bf_show, 'evidence')
    T = 2 * log(Q.bf_threshold);
    v(abs(v) < T) = 0;
end
out.dat(:, 1) = v;
if isa(out, 'statistic_image'), out.sig = v ~= 0; end

end


function m = brain_mask()

f = which('brainmask_canlab.nii');
if isempty(f), f = which('brainmask.nii'); end
if isempty(f)
    error('LaBGAScore_brain_figure:noBrainMask', 'brain_outline needs brainmask_canlab.nii (CanlabCore) on the path.');
end
m = fmri_mask_image(f, 'noverbose');

end


function out = apply_threshold(map, Q)

out = map;
if ~isempty(Q.threshold) && isa(map, 'statistic_image')
    out = threshold(map, Q.threshold{:}, 'noverbose');
end

end


function tf = is_thresholded(map)

tf = isa(map, 'statistic_image') && ~isempty(map.sig) && ~all(map.sig(:));

end


function v = shown_values(map)
% Values of the voxels that will actually be drawn.

if isa(map, 'region')
    v = [];
    for i = 1:numel(map)
        if ~isempty(map(i).Z), v = [v; double(map(i).Z(:))]; else, v = [v; double(map(i).val(:))]; end %#ok<AGROW>
    end
    v = v(isfinite(v) & v ~= 0);
    return
end

v = double(map.dat(:, 1));
if isa(map, 'statistic_image') && ~isempty(map.sig)
    v = v(logical(map.sig(:, 1)));              % .sig is stored as double 0/1
end
v = v(isfinite(v) & v ~= 0);

end


function v = shown_values_all(map)

v = double(map.dat(:));
v = v(isfinite(v) & v ~= 0);

end



%% =======================================================================
%  FIGURE, GRID AND SIZES
%  =======================================================================

function [fig, G] = open_figure(panels, P)

W = figure_width_mm(P.width);
n = numel(panels);

grid = P.grid;
if isempty(grid)
    if strcmp(P.layout, 'small_multiples')
        nc = P.ncols; if isempty(nc), nc = ceil(sqrt(n)); end
        grid = [ceil(n / nc) nc];
    else
        grid = [n 1];
    end
end
if prod(grid) < n
    error('LaBGAScore_brain_figure:grid', '''grid'' [%d %d] has fewer cells than the %d panels.', grid, n);
end

widths = P.widths; if isempty(widths), widths = ones(1, grid(2)); end
widths = widths / sum(widths);

% space reserved around the panels, in mm
has_title = ~isempty(P.title);
top = 1 + has_title * (P.fontsize * 0.5);
key_h = 0;
if needs_global_key(panels, P), key_h = 5 + P.fontsize * 0.9; end
side = 2;

% natural height of each row, from the panel types in it
row_h = zeros(1, grid(1));
for r = 1:grid(1)
    for c = 1:grid(2)
        i = (r - 1) * grid(2) + c;
        if i > n, continue, end
        cell_w = (W - 2 * side) * widths(c);
        row_h(r) = max(row_h(r), natural_height(panels{i}, cell_w, P));
    end
end
if ~isempty(P.heights)
    row_h = P.heights / sum(P.heights) * sum(row_h);
end

H = P.height;
if isempty(H), H = top + sum(row_h) + key_h + 1; end
scale = (H - top - key_h - 1) / sum(row_h);
row_h = row_h * scale;

fig = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2 2 W/10 H/10], ...
    'Visible', onoff(P.visible), 'InvertHardcopy', 'off', 'PaperPositionMode', 'auto', ...
    'Name', 'LaBGAScore_brain_figure');
set(fig, 'DefaultAxesFontName', P.fontname, 'DefaultTextFontName', P.fontname, ...
    'DefaultAxesFontSize', P.fontsize, 'DefaultTextFontSize', P.fontsize);

G = struct('W', W, 'H', H, 'rects', {cell(1, n)}, 'key_rect', [], 'grid', grid);
G.mm2x = @(mm) mm / W;
G.mm2y = @(mm) mm / H;

y_top = H - top;
for r = 1:grid(1)
    x = side;
    for c = 1:grid(2)
        i = (r - 1) * grid(2) + c;
        cw = (W - 2 * side) * widths(c);
        if i <= n
            G.rects{i} = [x / W, (y_top - row_h(r)) / H, cw / W, row_h(r) / H];
        end
        x = x + cw;
    end
    y_top = y_top - row_h(r);
end
if key_h > 0
    G.key_rect = [0.3, 1 / H, 0.4, (key_h - 1) / H];
end

end


function W = figure_width_mm(w)

if isnumeric(w), W = w; return, end
switch lower(char(w))
    case 'single',  W = 89;
    case 'onehalf', W = 120;
    case 'double',  W = 183;
    otherwise, error('LaBGAScore_brain_figure:width', '''width'' must be single, onehalf, double or a number in mm.');
end

end


function h = natural_height(pn, w, P)
% Height in mm a panel of width w needs to look right.

Q = pn.Q;
head = panel_head_mm(Q) + wants_labels(pn) * label_band_mm(Q);
foot = Q.slice_coords * (P.fontsize * 0.45 + 0.5);

switch pn.type
    case 'slices'
        rows = slice_rows(Q);
        h = head;
        for r = 1:numel(rows)
            [~, row_aspect] = row_geometry(rows{r});
            h = h + w / row_aspect + foot;
        end
    case 'surfaces'
        nv = numel(Q.views);
        h = head + w / nv * 0.72;
    case 'cutaway'
        h = head + min(w * 0.55, 70);
    case 'subcortex'
        nview = size(Q.view3d, 1); if nview == 0, nview = 2; end
        h = head + min(w / nview * 0.55, 40);
    case 'regioncenters'
        nr = numel(region_list(pn.map, Q));
        nr = min(max(nr, 1), Q.max_regions);
        nc = min(nr, 6);
        title_mm = 3 * (Q.fontsize - 1) * 0.3528 * 1.25 + 0.8;   % up to three title lines
        h = head + ceil(nr / nc) * (w / nc + title_mm);
    case 'zoom'
        nmain = zoom_count(pn.map, Q);
        if nmain == 1
            h = max(head + w * 0.5 + 4, 30);
        else
            % a slice below, an inset above; both keep their aspect, so a
            % wider column only adds white space
            h = head + min(w / nmain, 55) * 1.7 + 6;
        end
    case 'topography'
        h = head + min(w * 0.75, 90);
    case 'small_multiple'
        kind = Q.panel{1};
        switch kind
            case 'axial',    a = 1.2;
            case 'sagittal', a = 0.83;
            case 'coronal',  a = 1.0;
            otherwise,       a = 0.72;
        end
        h = head + w * a;
    otherwise
        h = head + w * 0.4;
end

end


function tf = needs_global_key(panels, P)

tf = false;
if strcmpi(P.colorbar, 'none'), return, end
for i = 1:numel(panels)
    if ~panels{i}.own_map && ~strcmp(panels{i}.type, 'topography'), tf = true; return, end
end

end


function ax = make_overlay_axes(fig)
% One invisible full-figure axes in normalized units, for letters, labels,
% leader lines and legends.

ax = axes('Parent', fig, 'Position', [0 0 1 1], 'Visible', 'off', 'XLim', [0 1], 'YLim', [0 1], ...
    'HitTest', 'off', 'Tag', 'LaBGAScore_overlay');
hold(ax, 'on');

end


function s = onoff(tf)

if tf, s = 'on'; else, s = 'off'; end

end



%% =======================================================================
%  DRAWING PANELS
%  =======================================================================

function [pn, o] = draw_panel(pn, o, fig, rect, ~)

Q = pn.Q;
pn.rect = rect;
H = fig_mm(fig);
head = panel_head_mm(Q) / H(2);
inner = [rect(1), rect(2), rect(3), rect(4) - head];

% panel letter and title
idx = panel_index(fig);
if ~(islogical(Q.panel_letters) && ~Q.panel_letters) && ~isequal(Q.panel_letters, false)
    letter = char(double(char(Q.panel_letters)) + idx - 1);
    annotation(fig, 'textbox', [rect(1), rect(2) + rect(4) - head, 0.05, head], 'String', letter, ...
        'EdgeColor', 'none', 'FontWeight', 'bold', 'FontSize', Q.fontsize + 2, 'FontName', Q.fontname, ...
        'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', 'Margin', 0);
end
if ~isempty(Q.panel_title)
    annotation(fig, 'textbox', [rect(1), rect(2) + rect(4) - head, rect(3), head], 'String', Q.panel_title, ...
        'EdgeColor', 'none', 'FontSize', Q.fontsize, 'FontName', Q.fontname, 'Interpreter', 'none', ...
        'VerticalAlignment', 'top', 'HorizontalAlignment', 'center', 'Margin', 0);
end
set_panel_index(fig, idx + 1);

% room for a panel's own colorbar when its map differs from the global one
if pn.own_map && ~strcmp(pn.type, 'topography') && ~strcmpi(Q.colorbar, 'none')
    kw = 14 / H(1);
    pn.key_rect = [inner(1) + inner(3) - kw, inner(2), kw, inner(4)];
    inner(3) = inner(3) - kw;
end

% a band above the content for labels, so they never sit on the brain
pn.label_band = [];
if wants_labels(pn)
    band = label_band_mm(Q) / H(2);
    pn.label_band = [inner(1), inner(2) + inner(4) - band, inner(3), band];
    inner(4) = inner(4) - band;
end

switch pn.type
    case 'slices',         [pn, o] = draw_slices(pn, o, fig, inner);
    case 'surfaces',       [pn, o] = draw_surfaces(pn, o, fig, inner);
    case 'cutaway',        [pn, o] = draw_cutaway(pn, o, fig, inner);
    case 'subcortex',      [pn, o] = draw_subcortex(pn, o, fig, inner);
    case 'regioncenters',  [pn, o] = draw_regioncenters(pn, o, fig, inner);
    case 'zoom',           [pn, o] = draw_zoom(pn, o, fig, inner);
    case 'topography',     pn = draw_topography(pn, fig, inner);
    case 'small_multiple', [pn, o] = draw_small_multiple(pn, o, fig, inner);
    otherwise
        error('LaBGAScore_brain_figure:badLayout', 'Unknown layout or panel type ''%s''.', pn.type);
end

end


function tf = wants_labels(pn)

tf = ~isequal(pn.Q.labels, false) && ~any(strcmp(pn.type, {'regioncenters', 'topography'}));

end


function h = label_band_mm(Q)
% Two lines of labels, alternating, so neighbouring names do not collide.

h = 2 * (Q.fontsize - 1) * 0.3528 * 1.4 + 1.5;

end


function h = panel_head_mm(Q)
% Space above a panel's content: room for its letter / title, if it has one.

has_letter = ~(islogical(Q.panel_letters) && ~Q.panel_letters) && ~isequal(Q.panel_letters, false);
h = 0.5 + (has_letter || ~isempty(Q.panel_title)) * (Q.fontsize * 0.5 + 0.5);

end


function idx = panel_index(fig)

idx = getappdata(fig, 'lbf_panel_index');
if isempty(idx), idx = 1; end

end


function set_panel_index(fig, idx)

setappdata(fig, 'lbf_panel_index', idx);

end


function mm = fig_mm(fig)

u = get(fig, 'Units');
set(fig, 'Units', 'centimeters');
pos = get(fig, 'Position');
set(fig, 'Units', u);
mm = pos(3:4) * 10;

end


% ---------------------------------------------------------------- slices

function rows = slice_rows(Q)

rows = Q.slices;
if isempty(rows)
    % x = -2 as in CanlabCore's compact2: x = 0 falls in the interhemispheric
    % fissure, where the template has no occipital cortex
    rows = {{{'sagittal', -2}, {'axial', -30:10:60}}};
end
for r = 1:numel(rows)
    if ischar(rows{r}{1}), rows{r} = rows(r); end       % single segment -> list of one
end

end


function [slices, aspect] = row_geometry(row)
% Slices of a row and the row's width/height ratio, from the MNI extent of
% each orientation (x -90..90, y -126..90, z -72..108 mm).

ext = struct('axial', [180 216], 'coronal', [180 180], 'sagittal', [216 180]);
slices = struct('orientation', {}, 'coord', {}, 'w', {});
for s = 1:numel(row)
    orient = normalise_orientation(row{s}{1});
    coords = row{s}{2};
    e = ext.(orient);
    for c = 1:numel(coords)
        slices(end+1) = struct('orientation', orient, 'coord', coords(c), 'w', e(1) / e(2)); %#ok<AGROW>
    end
end
aspect = sum([slices.w]) * 1.04;            % 4% gap between slices

end


function o = normalise_orientation(o)

o = lower(char(o));
if any(strcmp(o, {'saggital', 'sagg', 'sag'})), o = 'sagittal'; end
if strcmp(o, 'cor'), o = 'coronal'; end
if strcmp(o, 'ax'), o = 'axial'; end

end


function xyz = slice_xyz(orientation, coord)

switch orientation
    case 'axial',    xyz = [0 0 coord];
    case 'coronal',  xyz = [0 coord 0];
    case 'sagittal', xyz = [coord 0 0];
end

end


function o = add_montage(o, fig, orientation, xyz, ax)
% One slice in one existing axes. The underlay volume is read once per figure
% and handed back to montage ('volume_data'), instead of being re-read for
% every slice - which is most of the run time of a 12-slice row otherwise.

dat = getappdata(fig, 'lbf_underlay');
if isempty(dat)
    [o, dat] = montage(o, orientation, 'wh_slice', xyz, 'existing_axes', ax, 'noverbose');
    setappdata(fig, 'lbf_underlay', dat);
else
    o = montage(o, 'volume_data', dat, orientation, 'wh_slice', xyz, 'existing_axes', ax, 'noverbose');
end

end


function [pn, o] = draw_slices(pn, o, fig, inner)

Q = pn.Q;
rows = slice_rows(Q);
H = fig_mm(fig);
foot = Q.slice_coords * (Q.fontsize * 0.45 + 0.5) / H(2);

if Q.slice_locator
    loc_w = inner(4) / numel(rows) * H(2) / H(1) * 1.0;   % a square-ish sagittal
    loc_rect = [inner(1), inner(2), loc_w, inner(4)];
    inner(1) = inner(1) + loc_w * 1.05;
    inner(3) = inner(3) - loc_w * 1.05;
end

row_h = inner(4) / numel(rows);
for r = 1:numel(rows)
    [sl, aspect] = row_geometry(rows{r});
    y0 = inner(2) + inner(4) - r * row_h + foot;
    h = row_h - foot;
    % fit the row into the rect: limited by width or by height
    total_w = h * H(2) * aspect / H(1);
    if total_w > inner(3)
        h = h * inner(3) / total_w;
        total_w = inner(3);
    end
    x = inner(1) + (inner(3) - total_w) / 2;
    for s = 1:numel(sl)
        w = total_w * sl(s).w / sum([sl.w]) / 1.04;
        ax = axes('Parent', fig, 'Position', [x, y0, w, h]);
        o = add_montage(o, fig, sl(s).orientation, slice_xyz(sl(s).orientation, sl(s).coord), ax);
        pn.montages(end+1) = numel(o.montage);
        pn.axes(end+1) = ax;
        pn.slice_info(end+1) = struct('ax', ax, 'orientation', sl(s).orientation, 'coord', sl(s).coord);
        x = x + w * 1.04;
    end
end

if Q.slice_locator
    [pn, o] = draw_slice_locator(pn, o, fig, loc_rect);
end

end


function [pn, o] = draw_slice_locator(pn, o, fig, rect)
% Midline sagittal with a line at each axial / coronal slice (Sun et al. 2026
% Fig 5). It is a montage of o, but not in pn.montages, so it gets no blobs.

ax = axes('Parent', fig, 'Position', rect);
o = add_montage(o, fig, 'sagittal', [-2 0 0], ax);
hold(ax, 'on');
xl = [-105 72]; zl = [-75 85];
for k = 1:numel(pn.slice_info)
    s = pn.slice_info(k);
    switch s.orientation
        case 'axial',   plot(ax, xl, [s.coord s.coord], '-', 'Color', [0.1 0.3 0.8], 'LineWidth', 0.4);
        case 'coronal', plot(ax, [s.coord s.coord], zl, '-', 'Color', [0.1 0.3 0.8], 'LineWidth', 0.4);
    end
end

end


% ---------------------------------------------------------------- surfaces

function [pn, o] = draw_surfaces(pn, o, ~, inner)

Q = pn.Q;
nv = numel(Q.views);
w = inner(3) / nv;
for v = 1:nv
    [kw, orient, fslr] = surface_keyword(Q.surface, Q.views{v});
    pos = [inner(1) + (v - 1) * w, inner(2), w, inner(4)];
    n0 = numel(o.surface);
    o = surface(o, 'direction', kw, 'orientation', orient, 'axes', pos);
    new = n0 + 1:numel(o.surface);
    if fslr && Q.projection
        pn.surf_fslr = [pn.surf_fslr new];
    else
        pn.surf_mni = [pn.surf_mni new];
    end
    ax = surface_axes(o, new);
    pn.view_axes = [pn.view_axes ax];
end
same_scale(pn.view_axes);

end


function same_scale(axs)
% addbrain frames each view to fit its own mesh, so a lateral and a medial view
% of the same hemisphere come out at different sizes. One camera view angle -
% the widest any view needs - gives every view the same scale.

axs = axs(isgraphics(axs));
if numel(axs) < 2, return, end
ang = zeros(1, numel(axs));
for k = 1:numel(axs)
    set(axs(k), 'CameraViewAngleMode', 'auto');
    ang(k) = get(axs(k), 'CameraViewAngle');
end
set(axs, 'CameraViewAngle', max(ang));

end


function [kw, orient, fslr] = surface_keyword(kind, view_name)

view_name = lower(view_name);
hemi = 'left';
if contains(view_name, 'right'), hemi = 'right'; end
orient = 'lateral';
if contains(view_name, 'medial'), orient = 'medial'; end

switch lower(char(kind))
    case 'inflated',     kw = ['hcp inflated ' hemi]; fslr = true;
    case 'flat',         kw = ['flat ' hemi]; fslr = true; orient = 'lateral';
    case 'midthickness', kw = ['MNI152NLin2009cAsym midthickness ' hemi]; fslr = false;
    otherwise
        error('LaBGAScore_brain_figure:surface', '''surface'' must be inflated, flat or midthickness.');
end

end


function ax = surface_axes(o, idx)

ax = gobjects(0);
for k = idx
    s = o.surface{k};
    if isfield(s, 'axis_handles') && ~isempty(s.axis_handles)
        ax = [ax s.axis_handles(:)']; %#ok<AGROW>
    elseif isfield(s, 'object_handle') && ~isempty(s.object_handle)
        ax = [ax ancestor(s.object_handle(1), 'axes')]; %#ok<AGROW>
    end
end
ax = unique(ax);

end


% ---------------------------------------------------------------- 3-D

function [pn, o] = draw_cutaway(pn, o, fig, inner)

Q = pn.Q;
before = findobj(fig, 'Type', 'axes');
n0 = numel(o.surface);
o = surface(o, 'direction', Q.cutaway, 'axes', inner);
pn.surf_mni = n0 + 1:numel(o.surface);

% some keywords (coronal_slabs_*) draw into several axes of their own; tile
% them across the panel rather than leaving them where addbrain put them
created = setdiff(findobj(fig, 'Type', 'axes'), before);
created = created(arrayfun(@(a) ~isempty(findobj(a, 'Type', 'patch')), created));
if numel(created) > 1
    w = inner(3) / numel(created);
    [~, order] = sort(arrayfun(@(a) a.Position(1), created));
    created = created(order);
    for k = 1:numel(created)
        set(created(k), 'Position', [inner(1) + (k - 1) * w, inner(2), w, inner(4)]);
    end
end
pn.view_axes = created(:)';
if ~isempty(Q.view3d)
    for a = pn.view_axes, view(a, Q.view3d); end
end
% some keywords (the insula and accumbens slabs) camzoom for a full-window
% figure; give the camera back its automatic fit so the mesh stays inside its
% panel instead of spilling over the neighbours and the colour key
for a = pn.view_axes
    axis(a, 'image', 'off');
    set(a, 'CameraViewAngleMode', 'auto');
end
same_scale(pn.view_axes);

end


function [pn, o] = draw_subcortex(pn, o, fig, inner)

Q = pn.Q;
ax = axes('Parent', fig, 'Position', inner);
pn.view_axes = ax;

if isa(pn.map, 'atlas')
    % an atlas is drawn as its own parcels, in the indexed colours
    axes(ax);
    h = isosurface(pn.map);
    axis(ax, 'image', 'vis3d', 'off'); material(ax, 'dull');
    v = Q.view3d; if isempty(v), v = [210 20]; end
    view(ax, v); lightRestoreSingle(ax);
    pn.layers = {};
    pn.isosurfaces = h;
    return
end

% one axes per view; the scene is rebuilt in each, because a graphics object
% can only have one parent
views = Q.view3d; if isempty(views), views = [-120 20; 120 20]; end
delete(ax);
pn.view_axes = gobjects(0);
nv = size(views, 1);
mm = fig_mm(fig);
w = min(inner(3) / nv, inner(4) * 1.25 * mm(2) / mm(1));     % about square views...
x0 = inner(1) + (inner(3) - nv * w) / 2;                      % ...centred in the panel
for k = 1:nv
    ax = axes('Parent', fig, 'Position', [x0 + (k - 1) * w, inner(2), w, inner(4)]);
    for s = 1:numel(Q.structures)
        n0 = numel(o.surface);
        o = surface(o, 'direction', Q.structures{s}, 'axes', ax);
        pn.surf_mni = [pn.surf_mni n0 + 1:numel(o.surface)];
    end
    axis(ax, 'image', 'vis3d', 'off');
    view(ax, views(k, :));
    lightRestoreSingle(ax);
    material(ax, 'dull');
    pn.view_axes(end+1) = ax;
end
same_scale(pn.view_axes);

end


% ---------------------------------------------------------------- region close-ups

function [pn, o] = draw_regioncenters(pn, o, fig, inner)
% One axis per region, centred on it and zoomed to its extent. Built here
% rather than with region/montage, which opens its own figure and so cannot
% be placed in a composite.

Q = pn.Q;
r = region_list(pn.map, Q);
nr = numel(r);
if nr == 0
    error('LaBGAScore_brain_figure:noRegions', 'No suprathreshold regions to show.');
end

if nr > Q.max_regions
    peak = arrayfun(@(x) max(abs(double(x.Z(:)))), r);
    [~, order] = sort(peak, 'descend');
    r = r(sort(order(1:Q.max_regions)));
    lbf_warn('LaBGAScore_brain_figure:manyRegions', ...
        'Showing the %d strongest of %d regions (''max_regions''). Pass Inf to show all.', Q.max_regions, nr);
    nr = numel(r);
end

names = region_names(r, Q);
pn.regions = r;

nc = min(nr, 6);
nrow = ceil(nr / nc);
H = fig_mm(fig);
title_h = (3 * (Q.fontsize - 1) * 0.3528 * 1.25 + 0.8) / H(2);   % up to three title lines
cw = inner(3) / nc;
ch = inner(4) / nrow;
orient = normalise_orientation(Q.orientation);

for k = 1:nr
    rr = ceil(k / nc); cc = k - (rr - 1) * nc;
    pos = [inner(1) + (cc - 1) * cw, inner(2) + inner(4) - rr * ch, cw * 0.96, ch - title_h];
    ax = axes('Parent', fig, 'Position', pos);
    c = r(k).mm_center;
    o = add_montage(o, fig, orient, c(:)', ax);
    pn.montages(end+1) = numel(o.montage);
    pn.axes(end+1) = ax;
    [u, v] = plane_axes(orient);
    xyz = r(k).XYZmm;
    pad = 8;
    pn.limits(end+1) = struct('ax', ax, ...
        'xlim', [min(xyz(u, :)) - pad, max(xyz(u, :)) + pad], ...
        'ylim', [min(xyz(v, :)) - pad, max(xyz(v, :)) + pad]);
    % drawn on the overlay in annotate_panel, above the zoomed data box: axes
    % titles came out about twice their set size headless and collided
    max_chars = max(8, floor(cw * 0.96 * H(1) / (0.62 * (Q.fontsize - 1) * 0.3528)));
    pn.titles(end+1) = struct('ax', ax, 'lines', {[wrap_words(names{k}, max_chars); {sprintf('(%d, %d, %d)', round(c))}]});
end

end


function lines = wrap_words(str, n)
% Break a label into lines of at most n characters, at spaces where possible.

words = strsplit(str, ' ');
lines = {};
cur = '';
for i = 1:numel(words)
    if isempty(cur)
        cur = words{i};
    elseif numel(cur) + 1 + numel(words{i}) <= n
        cur = [cur ' ' words{i}]; %#ok<AGROW>
    else
        lines{end+1, 1} = cur; %#ok<AGROW>
        cur = words{i};
    end
end
lines{end+1, 1} = cur;
lines = lines(1:min(2, end));                          % two lines of name at most

end


function r = region_list(map, Q)

if isa(map, 'region') || isa(map, 'atlas')
    r = to_region(map);
else
    r = region(apply_threshold(map, Q), 'noverbose');
end
r = r(arrayfun(@(x) ~isempty(x.XYZmm), r));

end


function names = region_names(r, Q)
% Name each region by the atlas label at its PEAK voxel (or the nearest
% labelled voxel within 6 mm), from the label level 'label_field'. The finest
% level of canlab2024 ('labels') gives names like 'BStem_Shen_Pons_Rcd'; the
% coarser levels read better in a figure.

names = arrayfun(@(k) sprintf('Region %d', k), 1:numel(r), 'UniformOutput', false);
if isequal(Q.labels, false), return, end
try
    [axyz, aidx, labs] = atlas_lookup(Q.label_atlas, Q.label_field);
    for k = 1:numel(r)
        z = double(r(k).Z(:));
        [~, ip] = max(abs(z));
        p = r(k).XYZmm(:, ip)';
        [d, j] = min(sum((axyz - p) .^ 2, 2));
        if sqrt(d) <= 6 && aidx(j) >= 1 && aidx(j) <= numel(labs)
            names{k} = strrep(labs{aidx(j)}, '_', ' ');
        end
    end
catch e
    lbf_warn('LaBGAScore_brain_figure:labels', 'Could not label regions from %s: %s', Q.label_atlas, e.message);
end

end


function [xyz, idx, labs] = atlas_lookup(name, field)
% Voxel coordinates, parcel indices and labels of an atlas, cached across
% calls because load_atlas takes seconds.

persistent CACHE
if isempty(CACHE), CACHE = containers.Map(); end
key = [name '|' field];
if ~isKey(CACHE, key)
    atl = load_atlas(name);
    if ~isprop(atl, field) || isempty(atl.(field))
        error('Atlas %s has no label level ''%s''.', name, field);
    end
    [xyz, dat] = vox_xyz(atl);
    keep = dat(:, 1) > 0;
    CACHE(key) = struct('xyz', xyz(keep, :), 'idx', round(dat(keep, 1)), 'labs', {atl.(field)});
end
c = CACHE(key);
xyz = c.xyz; idx = c.idx; labs = c.labs;

end


function [xyz, dat] = vox_xyz(obj)
% mm coordinates of an image's voxels, in the row order of obj.dat. The data
% rows skip removed voxels, so xyzlist has to skip them too.

ijk = double(obj.volInfo.xyzlist);
rv = obj.removed_voxels;
if numel(rv) == size(ijk, 1), ijk = ijk(~rv, :); end
xyz = (obj.volInfo.mat * [ijk ones(size(ijk, 1), 1)]')';
xyz = xyz(:, 1:3);
dat = double(obj.dat);

end


function [u, v] = plane_axes(orientation)
% Which MNI axes (1 = x, 2 = y, 3 = z) a montage slice shows as its X and Y.

switch orientation
    case 'axial',    u = 1; v = 2;
    case 'coronal',  u = 1; v = 3;
    case 'sagittal', u = 2; v = 3;
end

end


% ---------------------------------------------------------------- zoom insets

function [pn, o] = draw_zoom(pn, o, fig, inner)
% Main slice(s) along the bottom, insets above, joined by lines (Kragel et
% al. 2023 Fig 1a). An inset is the boxed part of the slice magnified or,
% with 'zoom_inset', 'topography', its values as a 3-D surface.

Q = pn.Q;
Z = zoom_specs(pn.map, Q);
topo = strcmpi(Q.zoom_inset, 'topography');
if topo && (~isa(pn.map, 'image_vector') || isa(pn.map, 'atlas'))
    error('LaBGAScore_brain_figure:zoomTopoInput', 'topography insets need one image (statistic_image / fmri_data).');
end
if topo && any(arrayfun(@(z) ~isempty(z.roi_name), Z)) && isempty(which('tpaps'))
    error('LaBGAScore_brain_figure:noCurveFitting', ...
        'ROI topography insets fit a spline with tpaps (Curve Fitting Toolbox), which is not on the path.');
end

nmain = numel(Z);
H = fig_mm(fig);
key_w = topo * (0.07 * inner(3) + 0.02);              % the insets' own colorbar, right
cw = (inner(3) - key_w) / nmain;
pn.zoom_links = struct('main', {}, 'inset', {}, 'box', {}, 'topo', {});
pn.zoom_topo = topo;
insets = struct('ax', {}, 'T', {}, 'name', {});
rois = {};

for m = 1:nmain
    orient = normalise_orientation(Z(m).orientation);
    coord = Z(m).coord;
    boxes = Z(m).boxes;
    names = Z(m).names;
    if ~isempty(Z(m).roi), rois{end+1} = Z(m).roi; end %#ok<AGROW>

    foot = Q.slice_coords * (Q.fontsize * 0.45 + 0.5) / H(2);
    main_pos = [inner(1) + (m - 1) * cw, inner(2) + foot, cw, inner(4) * 0.6 - foot];
    ax_main = axes('Parent', fig, 'Position', main_pos);
    o = add_montage(o, fig, orient, slice_xyz(orient, coord), ax_main);
    pn.montages(end+1) = numel(o.montage);
    pn.axes(end+1) = ax_main;
    pn.slice_info(end+1) = struct('ax', ax_main, 'orientation', orient, 'coord', coord);

    nb = size(boxes, 1);
    iw = cw / nb;
    % insets left to right in the order their boxes appear on the page, so
    % the connecting lines of neighbouring insets cannot cross (a sagittal
    % slice runs right to left, so the box order in the call says nothing)
    px = zeros(nb, 1);
    for b = 1:nb
        c = data_to_fig(ax_main, [mean(boxes(b, 1:2)) mean(boxes(b, 3:4)) 0], fig);
        px(b) = c(1);
    end
    [~, slot] = sort(px); slot(slot) = 1:nb;
    % room for the names ABOVE the insets: inside, they covered the very
    % detail the inset is there to show
    th = ~isempty(names) * (Q.fontsize - 1) * 0.3528 * 1.6 / H(2);
    for b = 1:nb
        bx = boxes(b, :);
        if topo
            ih = inner(4) * 0.38 - th;
        else
            aspect = (bx(4) - bx(3)) / (bx(2) - bx(1));    % height / width in mm
            ih = min(inner(4) * 0.36 - th, iw * aspect * H(1) / H(2));
        end
        % a 3-D inset needs room on its left for the height ticks
        ml = 0.04 + 0.1 * topo;
        ipos = [inner(1) + (m - 1) * cw + (slot(b) - 1) * iw + iw * ml, inner(2) + inner(4) - th - ih, iw * (0.96 - ml), ih];
        ax_in = axes('Parent', fig, 'Position', ipos);
        name = '';
        if numel(names) >= b, name = names{b}; end
        if topo
            insets(end+1) = struct('ax', ax_in, 'T', zoom_topography(pn.map, Q, orient, coord, bx, Z(m).roi_name), ...
                'name', name); %#ok<AGROW>
        else
            o = add_montage(o, fig, orient, slice_xyz(orient, coord), ax_in);
            pn.montages(end+1) = numel(o.montage);
            pn.axes(end+1) = ax_in;
            pn.limits(end+1) = struct('ax', ax_in, 'xlim', bx(1:2), 'ylim', bx(3:4));
            if ~isempty(name)
                pn.titles(end+1) = struct('ax', ax_in, 'lines', {{name}});
            end
        end
        pn.zoom_links(end+1) = struct('main', ax_main, 'inset', ax_in, 'box', bx, 'topo', topo);
    end
end

% the ROIs outlined on the main slices
for k = 1:numel(rois)
    pn.layers{end+1} = layer(rois{k}, {'outline', 'color', [0 0 0], 'linewidth', Q.outline_width}, true);
end

if topo
    pn = draw_zoom_topographies(pn, insets, inner, key_w);
end

end


function pn = draw_zoom_topographies(pn, insets, inner, key_w)
% Every inset on ONE colour and height scale, so they can be compared, with
% one colorbar for them all at the right of the inset row.

Q = pn.Q;
allz = cell2mat(arrayfun(@(i) i.T.Z(:), insets(:), 'UniformOutput', false));
clim = Q.cmaprange;
if numel(clim) ~= 2
    m = max(abs(allz), [], 'omitnan');
    if isempty(m) || m == 0, m = 1; end
    clim = [-m m];
end
zl = [min(allz, [], 'omitnan') max(allz, [], 'omitnan')];
cm = get_colormap(Q.colormap);
floor_on = ~isempty(Q.floor) && Q.floor;
for i = 1:numel(insets)
    draw_topo_surface(insets(i).ax, insets(i).T, clim, zl, cm, Q, floor_on, true);
    if i > 1, set(insets(i).ax, 'ZTickLabel', {}); end          % one shared scale, read on the first
    if ~isempty(insets(i).name)
        title(insets(i).ax, insets(i).name, 'FontSize', Q.fontsize - 1, 'FontWeight', 'normal', 'Interpreter', 'none');
    end
end

lab = Q.zlabel;
if strcmp(lab, 'value') && ~isempty(pn.C) && isfield(pn.C, 'label') && ~isempty(pn.C.label)
    lab = pn.C.label;
end
if startsWith(lab, '\it'), lab = ['{' lab '}']; end      % italic for the symbol only
pos = get(insets(end).ax, 'Position');
cb = colorbar(insets(end).ax, 'eastoutside');
cb.FontSize = Q.fontsize - 1;
cb.LineWidth = 0.4;
cb.Label.String = [lab ' (topographies)'];
cb.Label.FontSize = Q.fontsize - 1;
set(insets(end).ax, 'Position', pos);
cb.Position = [inner(1) + inner(3) - key_w + 0.01, pos(2) + pos(4) * 0.15, 0.01, pos(4) * 0.7];
pn.topography = struct('clim', clim, 'zlim', zl, 'n_insets', numel(insets));

end


function T = zoom_topography(map, Q, orient, coord, bx, roi_name)
% Boxes and clusters: the raw values of the boxed plane. ROIs: the spline
% fit, over the slice's own two axes and unfolded, so a bilateral ROI shows
% its two sides side by side, as on the slice below it.

if isempty(roi_name)
    Qs = Q; Qs.plane = {orient, coord}; Qs.box = bx;
    T = topography_slice(map, Qs);
else
    [u, v] = plane_axes(orient);
    names = {'x', 'y', 'z'};
    Qr = Q; Qr.roi = {roi_name}; Qr.projection_axes = names([u v]); Qr.hemisphere = 'both';
    T = topography_roi(map, Qr);
end

end


function Z = zoom_specs(map, Q)
% The 'zoom' option as main slices: given boxes, the strongest clusters of
% the map, or named ROIs.

Zin = Q.zoom;
if isempty(Zin)
    error('LaBGAScore_brain_figure:noZoom', 'layout ''zoom'' needs ''zoom''.');
end
if isstring(Zin) && ~isscalar(Zin), Zin = cellstr(Zin); end
Z = struct('orientation', {}, 'coord', {}, 'boxes', {}, 'names', {}, 'roi', {}, 'roi_name', {});

if ischarlike_(Zin)
    if ~strcmpi(Zin, 'clusters')
        error('LaBGAScore_brain_figure:zoomSpec', '''zoom'' as text must be ''clusters''.');
    end
    r = region_list(map, Q);
    if isempty(r)
        error('LaBGAScore_brain_figure:noClusters', ...
            '''zoom'', ''clusters'' found no clusters; pass a thresholded map or ''threshold''.');
    end
    peak = arrayfun(@(x) max(abs(double(x.Z(:)))), r);
    [~, ord] = sort(peak, 'descend');
    r = r(ord(1:min(Q.max_zoom, numel(r))));
    names = region_names(r, setfield(Q, 'labels', 'auto')); %#ok<SFLD>
    for k = 1:numel(r)
        [~, ip] = max(abs(double(r(k).Z(:))));
        Z(end+1) = zoom_spec(Q.orientation, r(k).XYZmm(:, ip)', r(k).XYZmm', names{k}, [], ''); %#ok<AGROW>
    end

elseif iscellstr(Zin) && ~any(strcmpi(Zin{1}, {'axial', 'sagittal', 'coronal'}))
    for k = 1:numel(Zin)
        roi = topography_roi_mask(setfield(Q, 'roi', Zin(k)), map); %#ok<SFLD>
        [xyz, dat] = vox_xyz(roi);
        xyz = xyz(dat(:, 1) ~= 0, :);
        Z(end+1) = zoom_spec(Q.orientation, mean(xyz, 1), xyz, strrep(Zin{k}, '_', ' '), roi, Zin{k}); %#ok<AGROW>
    end

else
    if ischar(Zin{1}), Zin = {Zin}; end
    for m = 1:numel(Zin)
        sp = Zin{m};
        nm = {};
        if numel(sp) >= 4, nm = sp{4}; end
        Z(end+1) = struct('orientation', sp{1}, 'coord', sp{2}, 'boxes', sp{3}, 'names', {nm}, ...
            'roi', [], 'roi_name', ''); %#ok<AGROW>
    end
end

end


function s = zoom_spec(orientation, centre, xyz, name, roi, roi_name)
% One main slice through centre, boxed 5 mm beyond xyz, at least 20 mm a side.

orient = normalise_orientation(orientation);
[u, v] = plane_axes(orient);
w = setdiff(1:3, [u v]);
lo = min(xyz(:, [u v]), [], 1) - 5;
hi = max(xyz(:, [u v]), [], 1) + 5;
c = (lo + hi) / 2;
half = max((hi - lo) / 2, 10);
box = [c(1) - half(1), c(1) + half(1), c(2) - half(2), c(2) + half(2)];
s = struct('orientation', orient, 'coord', round(centre(w)), 'boxes', box, 'names', {{name}}, ...
    'roi', roi, 'roi_name', roi_name);

end


function n = zoom_count(map, Q)
% Number of main slices, for the panel height.

if isstring(Q.zoom) && ~isscalar(Q.zoom), Q.zoom = cellstr(Q.zoom); end
if ischarlike_(Q.zoom)
    n = min(Q.max_zoom, max(1, numel(region_list(map, Q))));
elseif iscellstr(Q.zoom) && ~isempty(Q.zoom) && ~any(strcmpi(Q.zoom{1}, {'axial', 'sagittal', 'coronal'}))
    n = numel(Q.zoom);
elseif ~isempty(Q.zoom) && ischar(Q.zoom{1})
    n = 1;
else
    n = max(1, numel(Q.zoom));
end

end


function tf = ischarlike_(x)

tf = ischar(x) || isstring(x);

end


% ---------------------------------------------------------------- small multiples

function [pn, o] = draw_small_multiple(pn, o, fig, inner)

Q = pn.Q;
kind = lower(Q.panel{1});
if strcmp(kind, 'surface')
    view_name = 'left lateral';
    if numel(Q.panel) >= 2, view_name = Q.panel{2}; end
    [kw, orient, fslr] = surface_keyword(Q.surface, view_name);
    n0 = numel(o.surface);
    o = surface(o, 'direction', kw, 'orientation', orient, 'axes', inner);
    new = n0 + 1:numel(o.surface);
    if fslr && Q.projection, pn.surf_fslr = new; else, pn.surf_mni = new; end
    pn.view_axes = surface_axes(o, new);
else
    orient = normalise_orientation(kind);
    coord = 0;
    if numel(Q.panel) >= 2, coord = Q.panel{2}; end
    ax = axes('Parent', fig, 'Position', inner);
    o = add_montage(o, fig, orient, slice_xyz(orient, coord), ax);
    pn.montages(end+1) = numel(o.montage);
    pn.axes(end+1) = ax;
    pn.slice_info(end+1) = struct('ax', ax, 'orientation', orient, 'coord', coord);
end

end



%% =======================================================================
%  PAINTING
%  =======================================================================

function o = paint_panel(o, pn, ~)
% Each layer goes to this panel's montages and MNI-space surfaces in one
% addblobs call, and to its fsLR surfaces in a second call with projection.
% The two must be separate: render_on_surface switches projection OFF for a
% whole call as soon as ONE of its surfaces lacks a left/right tag, which
% every cutaway and subcortical patch does.

if isempty(pn.layers), return, end
proj = {'sourcespace', pn.Q.sourcespace, 'targetsurface', 'fsLR_32k'};

for i = 1:numel(pn.layers)
    L = pn.layers{i};
    if isempty_layer(L.src), continue, end
    src = to_region(L.src);
    mont = pn.montages;
    smni = pn.surf_mni;
    if L.montage_only, smni = []; end
    if ~isempty(mont) || ~isempty(smni)
        o = addblobs(o, src, L.args{:}, 'wh_montages', mont, 'wh_surfaces', smni, 'noverbose');
    end
    if ~isempty(pn.surf_fslr) && ~L.montage_only
        o = addblobs(o, src, L.args{:}, 'wh_montages', [], 'wh_surfaces', pn.surf_fslr, proj{:}, 'noverbose');
    end
end

end


function r = to_region(src)
% addblobs draws regions. An atlas must go through atlas2region - region()
% refuses atlases - which gives one region per parcel, valued by its index,
% exactly what 'indexmap' colours.

if isa(src, 'region')
    r = src;
elseif isa(src, 'atlas')
    r = atlas2region(src);
else
    r = region(src, 'noverbose');
end

end


function tf = isempty_layer(src)

if isa(src, 'region')
    tf = isempty(src) || all(arrayfun(@(x) isempty(x.XYZmm), src));
elseif isa(src, 'statistic_image') && ~isempty(src.sig)
    tf = ~any(src.sig(:) & src.dat(:) ~= 0);
else
    tf = ~any(src.dat(:) ~= 0);
end

end


function apply_axis_limits(pn)

for k = 1:numel(pn.limits)
    L = pn.limits(k);
    set(L.ax, 'XLim', sort(L.xlim), 'YLim', sort(L.ylim));
end

end



%% =======================================================================
%  ANNOTATION
%  =======================================================================

function annotate_panel(pn, overlay, fig, ~)

Q = pn.Q;

% slice coordinates and L / R
for k = 1:numel(pn.slice_info)
    s = pn.slice_info(k);
    if any(pn.axes == s.ax) && any(arrayfun(@(L) L.ax == s.ax, pn.limits))
        continue                                        % zoom insets carry no coordinate
    end
    % anchored to the slice's DATA box: 'axis image' letterboxes the slice
    % inside its axes, so the axes edges can be far from the brain
    if Q.slice_coords
        letter = struct('axial', 'z', 'coronal', 'y', 'sagittal', 'x');
        bb = data_box(s.ax, fig);
        text(overlay, (bb(1) + bb(3)) / 2, bb(2), sprintf('%s = %g', letter.(s.orientation), s.coord), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', Q.fontsize - 1);
    end
end
if Q.lr_labels
    first = find(arrayfun(@(s) any(strcmp(s.orientation, {'axial', 'coronal'})), pn.slice_info), 1);
    if ~isempty(first)
        bb = data_box(pn.slice_info(first).ax, fig);
        % below the slice in a zoom panel, where the insets fill the space above
        yy = bb(4); va = 'top';
        if strcmp(pn.type, 'zoom'), yy = bb(2); va = 'bottom'; end
        text(overlay, bb(1), yy, 'L', 'FontSize', Q.fontsize - 1, 'Color', [.4 .4 .4], ...
            'HorizontalAlignment', 'left', 'VerticalAlignment', va);
        text(overlay, bb(3), yy, 'R', 'FontSize', Q.fontsize - 1, 'Color', [.4 .4 .4], ...
            'HorizontalAlignment', 'right', 'VerticalAlignment', va);
    end
end

% region close-up titles, above each zoomed slice
for k = 1:numel(pn.titles)
    bb = data_box(pn.titles(k).ax, fig);
    text(overlay, (bb(1) + bb(3)) / 2, bb(4) + 0.004, pn.titles(k).lines, 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', 'FontSize', Q.fontsize - 1, 'Interpreter', 'none');
end

% zoom boxes and their connecting lines
if isfield(pn, 'zoom_links')
    for k = 1:numel(pn.zoom_links)
        z = pn.zoom_links(k);
        b = z.box;
        % the box on the overlay too, for the same reason as the markers
        bp = zeros(5, 2);
        bx = b([1 2 2 1 1]); by = b([3 3 4 4 3]);
        for j = 1:5, bp(j, :) = data_to_fig(z.main, [bx(j) by(j) 0], fig); end
        plot(overlay, bp(:, 1), bp(:, 2), 'k-', 'LineWidth', 0.6);
        if z.topo
            % a 3-D inset: join the box to the two outermost corners of its floor
            xl = get(z.inset, 'XLim'); yl = get(z.inset, 'YLim'); zf = min(get(z.inset, 'ZLim'));
            fc = zeros(4, 2);
            for j = 1:4
                fc(j, :) = data_to_fig(z.inset, [xl(1 + (j > 2)) yl(1 + mod(j, 2)) zf], fig);
            end
            [~, lo] = min(fc(:, 1)); [~, hi] = max(fc(:, 1));
            bc = [data_to_fig(z.main, [b(1) b(4) 0], fig); data_to_fig(z.main, [b(2) b(4) 0], fig)];
            ic = fc([lo hi], :);
            [~, ob] = sort(bc(:, 1));
            for j = 1:2
                plot(overlay, [bc(ob(j), 1) ic(j, 1)], [bc(ob(j), 2) ic(j, 2)], 'k-', 'LineWidth', 0.5);
            end
            continue
        end
        set(z.inset, 'Visible', 'on', 'XTick', [], 'YTick', [], 'Box', 'on', 'LineWidth', 0.6, 'XColor', 'k', 'YColor', 'k');
        % Join the box's upper corners to the inset's lower corners. Both are
        % projected from DATA coordinates: the axes Position is letterboxed by
        % 'axis image', and a sagittal slice runs right-to-left (XDir reverse),
        % so neither the Position nor the box's own corner order says which
        % corner is on the left of the page.
        bc = [data_to_fig(z.main, [b(1) b(4) 0], fig); data_to_fig(z.main, [b(2) b(4) 0], fig)];
        ic = [data_to_fig(z.inset, [b(1) b(3) 0], fig); data_to_fig(z.inset, [b(2) b(3) 0], fig)];
        [~, ob] = sort(bc(:, 1)); [~, oi] = sort(ic(:, 1));
        for j = 1:2
            plot(overlay, [bc(ob(j), 1) ic(oi(j), 1)], [bc(ob(j), 2) ic(oi(j), 2)], 'k-', 'LineWidth', 0.5);
        end
    end
end

% markers (small multiples)
if ~isempty(Q.markers)
    draw_markers(pn, overlay, fig, Q);
end

% labels with leader lines
if ~isequal(Q.labels, false) && ~strcmp(pn.type, 'regioncenters') && ~strcmp(pn.type, 'topography')
    draw_labels(pn, overlay, fig, Q);
end

end


function draw_markers(pn, overlay, fig, Q)

for m = 1:size(Q.markers, 1)
    xyz = Q.markers(m, :);
    for k = 1:numel(pn.slice_info)
        s = pn.slice_info(k);
        [u, v] = plane_axes(s.orientation);
        w = setdiff(1:3, [u v]);
        if abs(xyz(w) - s.coord) <= 6
            % on the overlay, which is always on top: drawn into the slice
            % axes, the ring can end up under CanlabCore's slice surfaces
            pt = data_to_fig(s.ax, [xyz(u) xyz(v) 0], fig);
            plot(overlay, pt(1), pt(2), 'o', 'MarkerSize', 6, 'LineWidth', 1, 'Color', Q.marker_color);
        end
    end
    for a = pn.view_axes
        [pt, visible] = surface_point(a, xyz, fig);
        if visible
            plot(overlay, pt(1), pt(2), 'o', 'MarkerSize', 6, 'LineWidth', 1, 'Color', Q.marker_color);
        end
    end
end

end


function draw_labels(pn, overlay, fig, Q)
% Region names in the band above the panel, each over its point and joined to
% it by a dotted leader line (Kwon et al. 2026 Fig 5a). Each name is drawn
% once per panel: at the nearest slice, or in the first view that shows it.

labels = label_list(pn, Q);
if isempty(labels) || isempty(pn.label_band), return, end

pts = zeros(0, 2); txt = {}; col = zeros(0, 3);
for i = 1:size(labels, 1)
    xyz = labels{i, 2};
    pt = [];
    % slices: the slice nearest the point, if within 5 mm
    best = []; bestd = 5;
    for k = 1:numel(pn.slice_info)
        s = pn.slice_info(k);
        if any(arrayfun(@(L) L.ax == s.ax, pn.limits)), continue, end
        [u, v] = plane_axes(s.orientation);
        w = setdiff(1:3, [u v]);
        d = abs(xyz(w) - s.coord);
        if d <= bestd, best = k; bestd = d; end
    end
    if ~isempty(best)
        s = pn.slice_info(best);
        [u, v] = plane_axes(s.orientation);
        pt = data_to_fig(s.ax, [xyz(u) xyz(v) 0], fig);
    end
    % surfaces and 3-D views
    for a = pn.view_axes
        if ~isempty(pt), break, end
        [p2, visible] = surface_point(a, xyz, fig);
        if visible, pt = p2; end
    end
    if ~isempty(pt)
        pts(end+1, :) = pt; txt{end+1} = labels{i, 1}; col(end+1, :) = labels{i, 3}; %#ok<AGROW>
    end
end
if isempty(pts), return, end

mm = fig_mm(fig);
fs = Q.fontsize - 1;
line_h = fs * 0.3528 * 1.4 / mm(2);
gap = 1.5 / mm(1);
B = pn.label_band;
row_y = B(2) + 0.75 / mm(2) + line_h * [0.5 1.5];       % lower, upper line
L = B(1); R = B(1) + B(3);

% measure every label first (floored at a generous estimate - 0.62 em per
% character - because Extent under-reports headless)
[~, ord] = sort(pts(:, 1));
n = numel(ord);
h = gobjects(n, 1); w = zeros(n, 1);
for i = 1:n
    j = ord(i);
    h(i) = text(overlay, 0, row_y(1), txt{j}, 'Color', col(j, :), 'FontSize', fs, 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', 'Interpreter', 'none');
    ext = get(h(i), 'Extent');
    w(i) = max(ext(3), numel(txt{j}) * 0.62 * fs * 0.3528 / mm(1));
end

% assign each label to the line where it needs the smallest shift
row = zeros(n, 1); right_edge = [-Inf -Inf];
for i = 1:n
    x0 = pts(ord(i), 1);
    shift = max(0, right_edge + gap + w(i) / 2 - x0);
    [~, row(i)] = min(shift);
    right_edge(row(i)) = max(x0, right_edge(row(i)) + gap + w(i) / 2) + w(i) / 2;
end

% per line: a forward pass pushes labels right past their neighbours, a
% backward pass pulls them back inside the right edge, pushing earlier ones
% left. (Clamping after a one-way shift put a label straight back on top of
% its neighbour at the right edge.)
x = zeros(n, 1);
for r = 1:2
    k = find(row == r);
    if isempty(k), continue, end
    xr = pts(ord(k), 1);
    xr(1) = max(xr(1), L + w(k(1)) / 2);
    for m = 2:numel(k)
        xr(m) = max(xr(m), xr(m - 1) + (w(k(m - 1)) + w(k(m))) / 2 + gap);
    end
    xr(end) = min(xr(end), R - w(k(end)) / 2);
    for m = numel(k) - 1:-1:1
        xr(m) = min(xr(m), xr(m + 1) - (w(k(m)) + w(k(m + 1))) / 2 - gap);
    end
    if xr(1) - w(k(1)) / 2 < L - 1e-6
        lbf_warn('LaBGAScore_brain_figure:crowdedLabels', ...
            'Labels do not fit in one panel''s band without overlap; lower ''max_labels'' or widen the figure.');
    end
    x(k) = xr;
end

for i = 1:n
    j = ord(i);
    set(h(i), 'Position', [x(i) row_y(row(i)) 0]);
    plot(overlay, [pts(j, 1) x(i)], [pts(j, 2) row_y(row(i)) - line_h * 0.45], ':', 'Color', [0.25 0.25 0.25], 'LineWidth', 0.5);
end
% names on top of the lines, on white: an upper-line leader passes through
% the lower line of names, and crossed out it read as two names overlapping
set(h, 'BackgroundColor', 'w', 'Margin', 0.1);
uistack(h, 'top');

end


function labels = label_list(pn, Q)
% n x 3 cell: name, [x y z], colour.

labels = cell(0, 3);
if iscell(Q.labels)
    for i = 1:size(Q.labels, 1)
        labels(end+1, :) = {Q.labels{i, 1}, Q.labels{i, 2}, [0 0 0]}; %#ok<AGROW>
    end
    return
end
if ~strcmpi(Q.labels, 'auto') || iscell(pn.map) || isa(pn.map, 'atlas') || isa(pn.map, 'fmri_mask_image')
    return
end

src = pn.map;
if strcmpi(Q.style, 'unthresholded_outline')
    [~, src] = split_significance(pn.map, Q);
end
r = region_list(src, Q);
if isempty(r), return, end
peak = arrayfun(@(x) max(abs(double(x.Z(:)))), r);
[~, ord] = sort(peak, 'descend');
r = r(ord(1:min(Q.max_labels, numel(r))));
names = region_names(r, setfield(Q, 'labels', 'auto')); %#ok<SFLD>
for k = 1:numel(r)
    z = double(r(k).Z(:));
    [~, ip] = max(abs(z));
    xyz = r(k).XYZmm(:, ip)';
    if z(ip) > 0, c = [0.75 0.25 0]; else, c = [0 0.3 0.75]; end
    labels(end+1, :) = {names{k}, xyz, c}; %#ok<AGROW>
end

end


function [pt, visible] = surface_point(ax, xyz, fig)
% Figure-normalized position of an MNI point on a surface view, and whether
% it faces the camera. fsLR meshes are not in MNI space: the point is first
% located on the MNI152NLin2009cAsym midthickness mesh and carried to the fsLR
% vertex the projection maps it to, using the same nearest-neighbour file the
% projection itself uses.

pt = [NaN NaN]; visible = false;
p = findobj(ax, 'Type', 'patch');
if isempty(p), return, end

best = Inf; vert = []; nrm = []; target = [];
for h = p(:)'
    tag = lower(get(h, 'Tag'));
    V = get(h, 'Vertices');
    if size(V, 1) == 32492 && (contains(tag, 'left') || contains(tag, 'right'))
        hemi = 'lh'; if contains(tag, 'right'), hemi = 'rh'; end
        [vi, d] = fslr_vertex(xyz, hemi);
    else
        [d, vi] = min(sum((V - xyz) .^ 2, 2));
        d = sqrt(d);
    end
    if d < best
        best = d; vert = V(vi, :); target = h;
        N = get(h, 'VertexNormals');
        if ~isempty(N), nrm = N(vi, :); end
    end
end
if isempty(target) || best > 6, return, end

campos = get(ax, 'CameraPosition');
if ~isempty(nrm)
    % VertexNormals point inward for many CANlab meshes; accept either sign
    % only when the vertex is on the silhouette side facing the camera
    facing = dot(nrm / norm(nrm), (campos - vert) / norm(campos - vert));
    visible = abs(facing) > 0.15 && sign(facing) == mode_normal_sign(target, campos);
else
    visible = true;
end
pt = data_to_fig(ax, vert, fig);

end


function s = mode_normal_sign(h, campos)
% Whether this mesh's normals point outward (+1) or inward (-1), judged from
% the vertices nearest the camera, which must face it.

V = get(h, 'Vertices'); N = get(h, 'VertexNormals');
if isempty(N), s = 1; return, end
d = sum((V - campos) .^ 2, 2);
[~, ord] = sort(d);
near = ord(1:min(200, numel(ord)));
f = sum(N(near, :) .* (campos - V(near, :)), 2);
s = sign(median(f));
if s == 0, s = 1; end

end


function [vi, d] = fslr_vertex(xyz, hemi)

persistent REG SRC
if isempty(REG)
    REG = importdata(which('resample_from_MNI152NLin2009cAsym_to_fsLR_32k_nearestneighbor.mat'));
    SRC.lh = load(which('MNI152NLin2009cAsym_midthickness_lh.mat'), 'vertices');
    SRC.rh = load(which('MNI152NLin2009cAsym_midthickness_rh.mat'), 'vertices');
end
srcv = SRC.(hemi).vertices;
map = REG.(['vertices_' hemi]);
pos = srcv(map, :);                              % MNI position each fsLR vertex samples
[d, vi] = min(sum((pos - xyz) .^ 2, 2));
d = sqrt(d);

end


function bb = data_box(ax, fig)
% [left bottom right top] in figure-normalized units of what a 2-D slice axes
% actually shows: its data limits, projected.

xl = get(ax, 'XLim'); yl = get(ax, 'YLim');
p1 = data_to_fig(ax, [xl(1) yl(1) 0], fig);
p2 = data_to_fig(ax, [xl(2) yl(2) 0], fig);
bb = [min(p1(1), p2(1)), min(p1(2), p2(2)), max(p1(1), p2(1)), max(p1(2), p2(2))];

end


function pt = data_to_fig(ax, xyz, fig)
% Data coordinates -> figure-normalized coordinates, through MATLAB's own
% camera transform, so it holds for 2-D slices and rotated 3-D views alike.

px = matlab.graphics.chart.internal.convertDataSpaceCoordsToViewerCoords(ax, xyz(:));
u = get(fig, 'Units');
set(fig, 'Units', 'pixels');
fp = get(fig, 'Position');
set(fig, 'Units', u);
pt = [px(1) / fp(3), px(2) / fp(4)];

end



%% =======================================================================
%  COLORBARS AND LEGENDS
%  =======================================================================

function info = draw_colour_keys(panels, overlay, G, fig, P, info)

drawn_global = false;
for i = 1:numel(panels)
    pn = panels{i};
    if isempty(pn.C) || strcmp(pn.type, 'topography'), continue, end
    leg = layer_legend(pn);
    if ~isempty(pn.key_rect)
        draw_key(fig, overlay, pn.C, leg, pn.key_rect, 'east', P);
    elseif ~drawn_global && ~isempty(G.key_rect)
        C = pn.C;
        C.label = strtrim([C.label ' ' key_panels(panels, P)]);
        draw_key(fig, overlay, C, leg, G.key_rect, 'south', P);
        drawn_global = true;
        if isfield(pn.C, 'bayes'), info.bayes = pn.C.bayes; end
    end
end

end


function s = key_panels(panels, P)
% '(A-C)' when the shared key does not cover every panel (a topography or a
% panel with its own key), so the reader knows which panels it belongs to.

s = '';
if numel(panels) == 1 && isfield(panels{1}, 'zoom_topo') && panels{1}.zoom_topo
    s = '(slices)';
    return
end
if numel(panels) < 2 || isequal(P.panel_letters, false), return, end
on = cellfun(@(pn) ~isempty(pn.C) && isempty(pn.key_rect) && ~strcmp(pn.type, 'topography'), panels);
if all(on), return, end
k = find(on);
letters = char(double(char(P.panel_letters)) + k - 1);
if numel(k) > 2 && all(diff(k) == 1)
    s = sprintf('(%s%s%s)', letters(1), char(8211), letters(end));
else
    s = ['(' strjoin(cellstr(letters(:)), ', ') ')'];
end

end


function leg = layer_legend(pn)
% Legend entries carried by the layers themselves (multi-threshold levels).

leg = cell(0, 2);
for i = 1:numel(pn.layers)
    if isfield(pn.layers{i}, 'legend'), leg = [leg; pn.layers{i}.legend]; end %#ok<AGROW>
end

end


function draw_key(fig, overlay, C, leg, rect, where, P)

if strcmpi(P.colorbar, 'none'), return, end

% swatch legend: categories, atlas regions, threshold levels, solid masks
if ~isempty(leg) || any(strcmp(C.mode, {'indexed', 'categorical', 'solid'}))
    if isempty(leg) && ~isempty(C.legend)
        leg = [{C.legend.name}' {C.legend.color}'];
    end
    if isempty(leg), return, end                       % e.g. an atlas of more than 24 regions
    if strcmp(where, 'south')
        rect = [0.03, rect(2), 0.94, rect(4)];          % swatches use the full width
    end
    draw_swatches(fig, overlay, leg, rect, P);
    return
end

% colorbar, drawn from the SAME canlab_colormap object the blobs were painted
% with, so it cannot disagree with them. A thresholded split map gets two
% ramps, negative and positive, with nothing for the sub-threshold gap.
if strcmp(C.mode, 'split') && C.range(3) > 0
    segments = {[C.range(1) C.range(2)], [C.range(3) C.range(4)]};
elseif strcmp(C.mode, 'split')
    segments = {[C.range(1) C.range(4)]};
elseif strcmp(C.mode, 'discrete')
    segments = {[C.range(1) - 0.5, C.range(2) + 0.5]};
else
    segments = {[C.range(1) C.range(end)]};
end

mm = fig_mm(fig);
fs = P.fontsize - 1;
nseg = numel(segments);
gap = 0.16;                                            % fraction of the key's length between segments,
                                                       % wide enough that the two inner tick labels do not meet

if strcmp(where, 'south')
    bar_h = 1.8 / mm(2);                               % 1.8 mm high
    top = rect(2) + rect(4);
    y_bar = top - (fs * 0.3528 * 1.7 + 0.8) / mm(2) - bar_h;   % room for the label above, subscripts included
    text(overlay, rect(1) + rect(3) / 2, top, C.label, 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'top', 'FontSize', fs, 'Interpreter', 'tex');
    seg_w = (rect(3) * (1 - gap * (nseg - 1))) / nseg;
    for k = 1:nseg
        lim = segments{k};
        x0 = rect(1) + (k - 1) * (seg_w + gap * rect(3));
        ax = axes('Parent', fig, 'Position', [x0, y_bar, seg_w, bar_h]);
        draw_ramp(ax, C, lim, 'h');
        [t, l] = segment_ticks(C, lim);
        set(ax, 'XTick', t, 'XTickLabel', l, 'FontSize', fs);
    end
else
    bar_w = 1.8 / mm(1);
    seg_h = (rect(4) * 0.7 * (1 - gap * (nseg - 1))) / nseg;
    x0 = rect(1) + rect(3) * 0.25;
    for k = 1:nseg
        lim = segments{k};
        y0 = rect(2) + rect(4) * 0.15 + (k - 1) * (seg_h + gap * rect(4) * 0.7);
        ax = axes('Parent', fig, 'Position', [x0, y0, bar_w, seg_h]);
        draw_ramp(ax, C, lim, 'v');
        [t, l] = segment_ticks(C, lim);
        set(ax, 'YTick', t, 'YTickLabel', l, 'FontSize', fs, 'YAxisLocation', 'right');
    end
    text(overlay, x0 + bar_w / 2, rect(2) + rect(4) * 0.88, C.label, 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', 'FontSize', fs, 'Interpreter', 'tex');
end

end


function draw_ramp(ax, C, lim, orientation)

n = 256;
if strcmp(C.mode, 'discrete')
    vals = (ceil(lim(1)):floor(lim(2)))';
else
    vals = linspace(lim(1), lim(2), n)';
end
rgb = C.cmap.map(vals);
rgb(any(isnan(rgb), 2), :) = 1;
if strcmp(orientation, 'h')
    image(ax, 'XData', [vals(1) vals(end)], 'YData', [0 1], 'CData', reshape(rgb, 1, [], 3));
    set(ax, 'XLim', lim, 'YLim', [0 1], 'YTick', []);
else
    image(ax, 'XData', [0 1], 'YData', [vals(1) vals(end)], 'CData', reshape(rgb, [], 1, 3));
    set(ax, 'YLim', lim, 'XLim', [0 1], 'XTick', [], 'YDir', 'normal');
end
% outline drawn as a rectangle: 'Box' on would mirror the ticks onto the far
% edge, where they poke into the label above the bar
set(ax, 'Box', 'off', 'LineWidth', 0.4, 'TickDir', 'out', 'TickLength', [0.02 0.02], 'Layer', 'top');
if strcmp(orientation, 'h')
    rectangle(ax, 'Position', [lim(1) 0 diff(lim) 1], 'EdgeColor', 'k', 'LineWidth', 0.4);
    set(ax, 'YColor', 'none');
else
    rectangle(ax, 'Position', [0 lim(1) 1 diff(lim)], 'EdgeColor', 'k', 'LineWidth', 0.4);
    set(ax, 'XColor', 'none');
end

end


function [t, l] = segment_ticks(C, lim)
% The colour key's ticks that fall on this segment.

keep = C.ticks >= min(lim) - 1e-9 & C.ticks <= max(lim) + 1e-9;
t = C.ticks(keep);
l = C.ticklabels(keep);
[t, ord] = unique(t);
l = l(ord);

end


function draw_swatches(fig, overlay, leg, rect, P)
% Colour swatch + name per entry, in as many equal columns as the longest
% name allows, filled row by row.

n = size(leg, 1);
mm = fig_mm(fig);
fs = P.fontsize - 1;
char_mm = 0.55 * fs * 0.3528;                          % average glyph width at fs pt
sw_mm = 2.2;
col_mm = sw_mm + 1 + max(cellfun(@numel, leg(:, 1))) * char_mm + 3;
ncol = max(1, min(n, floor(rect(3) * mm(1) / col_mm)));
nrow = ceil(n / ncol);
cw = col_mm / mm(1);
x0 = rect(1) + (rect(3) - ncol * cw) / 2;              % centred
rh = min(rect(4) / nrow, (fs * 0.3528 * 1.6) / mm(2));
y_top = rect(2) + rect(4) - rh / 2;
for k = 1:n
    rr = ceil(k / ncol); cc = k - (rr - 1) * ncol;
    x = x0 + (cc - 1) * cw;
    y = y_top - (rr - 1) * rh;
    rectangle(overlay, 'Position', [x, y - rh * 0.3, sw_mm / mm(1), rh * 0.6], 'FaceColor', leg{k, 2}, 'EdgeColor', 'none');
    text(overlay, x + (sw_mm + 1) / mm(1), y, leg{k, 1}, 'FontSize', fs, 'VerticalAlignment', 'middle', ...
        'Interpreter', 'none');
end

end



%% =======================================================================
%  TOPOGRAPHY
%  =======================================================================

function pn = draw_topography(pn, fig, inner)
% Cluster surface plot. 'roi': voxel values over a two-axis projection of an
% ROI, fitted with a thin-plate smoothing spline (Kragel et al. 2023 Fig 1c).
% 'slice': the raw values of one slice in a box (figures/cluster_surface_plots.m).

Q = pn.Q;
map = pn.map;
if iscell(map) || ~isa(map, 'image_vector')
    error('LaBGAScore_brain_figure:topoInput', 'topography needs one image (statistic_image / fmri_data).');
end

% With the ROI inset, a right-hand column holds the inset above the
% colorbar; overlaid on the plot, the inset covered the 3-D axes box.
has_inset = Q.inset && strcmpi(Q.topo_mode, 'roi');
key_w = 0.09 * inner(3) + 0.02;
if has_inset
    main_pos = [inner(1) + 0.08 * inner(3), inner(2) + 0.14 * inner(4), inner(3) * 0.62 - 0.02, inner(4) * 0.8];
    inset_pos = [inner(1) + 0.72 * inner(3), inner(2) + 0.52 * inner(4), inner(3) * 0.27, inner(4) * 0.44];
    cb_pos = [inner(1) + 0.8 * inner(3), inner(2) + 0.12 * inner(4), 0.012, inner(4) * 0.32];
else
    main_pos = [inner(1) + 0.08 * inner(3), inner(2) + 0.14 * inner(4), inner(3) * 0.74 - key_w, inner(4) * 0.8];
    cb_pos = [inner(1) + inner(3) - key_w, inner(2) + 0.3 * inner(4), 0.012, inner(4) * 0.45];
end
ax = axes('Parent', fig, 'Position', main_pos);

switch lower(Q.topo_mode)
    case 'roi',   T = topography_roi(map, Q);
    case 'slice', T = topography_slice(map, Q);
    otherwise, error('LaBGAScore_brain_figure:topoMode', '''topo_mode'' must be roi or slice.');
end

cm = get_colormap(Q.colormap);
clim = Q.cmaprange;
if isempty(clim)
    m = max(abs(T.Z(:)), [], 'omitnan');
    if isempty(m) || m == 0, m = 1; end
    clim = [-m m];
end

zl = [min(T.Z(:), [], 'omitnan') max(T.Z(:), [], 'omitnan')];
draw_topo_surface(ax, T, clim, zl, cm, Q, isempty(Q.floor) || Q.floor, false);
if ~isempty(T.title)
    title(ax, T.title, 'FontSize', Q.fontsize, 'FontWeight', 'normal', 'Interpreter', 'none');
end

% colorbar, from the surface's own colormap
cb = colorbar(ax, 'eastoutside');
cb.FontSize = Q.fontsize - 1;
cb.LineWidth = 0.4;
cb.Label.String = Q.zlabel;
cb.Label.FontSize = Q.fontsize - 1;
set(ax, 'Position', main_pos);
cb.Position = cb_pos;

if has_inset
    try
        draw_topography_inset(fig, inset_pos, T, map, cm, clim, Q);
    catch e
        lbf_warn('LaBGAScore_brain_figure:topoInset', 'ROI inset not drawn: %s', e.message);
    end
end

pn.axes = ax;
pn.topography = rmfield(T, intersect(fieldnames(T), {'Z', 'u', 'v', 'roi', 'spline'}));
pn.topography.clim = clim;

end


function T = topography_roi(map, Q)

[roi, level] = topography_roi_mask(Q, map);
[xyz, val] = roi_voxels(map, roi);

switch lower(Q.hemisphere)
    case 'left',  keep = xyz(:, 1) < 0;
    case 'right', keep = xyz(:, 1) > 0;
    otherwise,    keep = true(size(xyz, 1), 1);
end
xyz = xyz(keep, :); val = val(keep);
if numel(val) < 10
    error('LaBGAScore_brain_figure:topoFewVoxels', 'Only %d ROI voxels with non-zero values; need at least 10.', numel(val));
end

u = axis_values(xyz, Q.projection_axes{1});
v = axis_values(xyz, Q.projection_axes{2});

% folding the hemispheres (absx) or collapsing the third axis puts several
% voxels on one site; tpaps needs distinct sites, so average them first
[sites, ~, ic] = unique([u v], 'rows');
z = accumarray(ic, val, [], @mean);

if isempty(Q.smoothing)
    [st, p_used] = tpaps(sites', z');
else
    [st, p_used] = tpaps(sites', z', Q.smoothing);
end

n = 60;
ug = linspace(min(sites(:, 1)), max(sites(:, 1)), n);
vg = linspace(min(sites(:, 2)), max(sites(:, 2)), n);
Z = fnval(st, {ug, vg});

if isfinite(Q.mask_distance)
    [UU, VV] = ndgrid(ug, vg);
    d = sqrt(min((UU(:) - sites(:, 1)') .^ 2 + (VV(:) - sites(:, 2)') .^ 2, [], 2));
    Z(reshape(d, size(Z)) > Q.mask_distance) = NaN;
end

% descriptive linear gradients along each axis (not inference: voxels are not
% independent)
X = [ones(numel(val), 1) u v];
b = X \ val;

T = struct('mode', 'roi', 'u', ug, 'v', vg, 'Z', Z, 'axes', {Q.projection_axes}, 'spline', st, ...
    'n_voxels', numel(val), 'n_sites', size(sites, 1), 'smoothing_p', p_used, ...
    'slope_per_mm', b(2:3)', 'roi', roi, 'roi_label_level', level, 'hemisphere', lower(Q.hemisphere), ...
    'title', roi_title(Q));

end


function [roi, level] = topography_roi_mask(Q, map)

r = Q.roi;
level = '';
if isempty(r)
    error('LaBGAScore_brain_figure:noRoi', 'topography ''roi'' mode needs ''roi''.');
end
if ischar(r) || isstring(r), r = cellstr(r); end
if iscell(r)
    % labels live at different levels of an atlas (in canlab2024, 'VStriatum'
    % first appears at labels_3), so search the levels finest-first and use the
    % first at which EVERY requested label matches something
    atl = load_atlas(Q.roi_atlas);
    levels = {'labels', 'labels_2', 'labels_3', 'labels_4', 'labels_5'};
    roi = [];
    for f = levels
        if ~isprop(atl, f{1}) || isempty(atl.(f{1})), continue, end
        if all(cellfun(@(x) any(contains(atl.(f{1}), x)), r))
            roi = select_atlas_subset(atl, r, f{1}, 'flatten');
            level = f{1};
            break
        end
    end
    if isempty(roi) || num_regions(roi) == 0
        error('LaBGAScore_brain_figure:roiEmpty', 'No label level of %s matches all of: %s.', ...
            Q.roi_atlas, strjoin(r, ', '));
    end
elseif isa(r, 'region')
    roi = region2imagevec(r);
else
    roi = r;
end
roi = resample_space(roi, map, 'nearest');

end


function [xyz, val] = roi_voxels(map, roi)

mvol = unthresholded_volume(map);
rvol = reconstruct_image(roi);
idx = find(rvol(:) ~= 0 & mvol(:) ~= 0 & isfinite(mvol(:)));
[i, j, k] = ind2sub(size(mvol), idx);
xyz = (map.volInfo.mat * [i j k ones(numel(i), 1)]')';
xyz = xyz(:, 1:3);
val = double(mvol(idx));

end


function vol = unthresholded_volume(map)
% A topography shows the values AROUND a result, so it reads every voxel:
% reconstruct_image zeroes whatever a statistic_image's .sig excludes.

if isa(map, 'statistic_image') && ~isempty(map.sig)
    map.sig = true(size(map.sig));
end
vol = reconstruct_image(map);
vol = vol(:, :, :, 1);

end


function a = axis_values(xyz, name)

switch lower(name)
    case 'x',    a = xyz(:, 1);
    case 'absx', a = abs(xyz(:, 1));
    case 'y',    a = xyz(:, 2);
    case 'z',    a = xyz(:, 3);
    otherwise, error('LaBGAScore_brain_figure:projAxis', 'projection axes must be x, absx, y or z.');
end

end


function s = axis_label(name)

switch lower(name)
    case 'absx', s = '|x| (mm)   medial \leftrightarrow lateral';
    case 'x',    s = 'x (mm)   left \leftrightarrow right';
    case 'y',    s = 'y (mm)   posterior \leftrightarrow anterior';
    case 'z',    s = 'z (mm)   inferior \leftrightarrow superior';
    otherwise,   s = name;
end

end


function s = roi_title(Q)
% ROI name plus the hemispheres the fitted surface covers: with |x| both are
% folded onto one surface, which the two volumes of the inset do not show.

s = '';
if iscell(Q.roi) || ischar(Q.roi) || isstring(Q.roi)
    s = strjoin(cellstr(Q.roi), ', ');
    s = strrep(s, '_', ' ');
end
switch lower(Q.hemisphere)
    case 'left',  hemi = '(L)';
    case 'right', hemi = '(R)';
    otherwise
        if any(strcmpi(Q.projection_axes, 'absx')), hemi = '(L + R, folded at midline)';
        elseif any(strcmpi(Q.projection_axes, 'x')), hemi = '(L + R)';
        else, hemi = '(L + R, collapsed over x)';
        end
end
s = strtrim([s ' ' hemi]);
if ~isempty(Q.panel_title), s = ''; end

end


function T = topography_slice(map, Q)

if numel(Q.plane) < 2 || numel(Q.box) ~= 4
    error('LaBGAScore_brain_figure:topoSlice', '''slice'' mode needs ''plane'', {orientation, coord} and ''box'', [u1 u2 v1 v2].');
end
orient = normalise_orientation(Q.plane{1});
coord = Q.plane{2};

M = map.volInfo.mat;
if any(abs(M(1:3, 1:3) - diag(diag(M(1:3, 1:3)))) > 1e-6, 'all')
    error('LaBGAScore_brain_figure:oblique', 'slice topography needs an image without rotation in its affine.');
end
vol = unthresholded_volume(map);
dims = size(vol);
mm = cell(1, 3);
for d = 1:3
    mm{d} = M(d, d) * (1:dims(d)) + M(d, 4);
end

[u, v] = plane_axes(orient);
w = setdiff(1:3, [u v]);
[~, sl] = min(abs(mm{w} - coord));
switch w
    case 1, S = squeeze(vol(sl, :, :));
    case 2, S = squeeze(vol(:, sl, :));
    case 3, S = squeeze(vol(:, :, sl));
end
uu = mm{u}; vv = mm{v};
b = Q.box;
iu = uu >= min(b(1:2)) & uu <= max(b(1:2));
iv = vv >= min(b(3:4)) & vv <= max(b(3:4));
Z = double(S(iu, iv));
Z(Z == 0) = NaN;
[uu, ou] = sort(uu(iu)); [vv, ov] = sort(vv(iv));
Z = Z(ou, ov);

names = {'x', 'y', 'z'};
% which plane and which part of it: without this the plot has no anatomy
ttl = sprintf('%s slice %s = %s mm (%s %s to %s, %s %s to %s)', orient, names{w}, num_label(round(mm{w}(sl), 1)), ...
    names{u}, num_label(min(b(1:2))), num_label(max(b(1:2))), names{v}, num_label(min(b(3:4))), num_label(max(b(3:4))));
if ~isempty(Q.panel_title), ttl = ''; end
T = struct('mode', 'slice', 'u', uu, 'v', vv, 'Z', Z, 'axes', {names([u v])}, ...
    'n_voxels', nnz(isfinite(Z)), 'plane', {{orient, mm{w}(sl)}}, 'box', b, 'roi', [], 'title', ttl);

end


function draw_topo_surface(ax, T, clim, zl, cm, Q, floor_on, compact)
% The fitted or raw surface, on given colour (clim) and height (zl) limits so
% that several topographies can share one scale. compact: an inset - short
% axis labels, no height label.

hold(ax, 'on');
surf(ax, T.u, T.v, T.Z', 'EdgeColor', [0.15 0.15 0.15], 'EdgeAlpha', 0.25, 'FaceColor', 'interp', 'LineWidth', 0.2);
colormap(ax, cm);
caxis(ax, clim);

if diff(zl) == 0, zl = zl + [-1 1]; end
zfloor = zl(1) - 0.35 * diff(zl) * floor_on;
if floor_on
    draw_topography_floor(ax, T, zfloor, cm, clim);
end

set(ax, 'ZLim', [zfloor zl(2) + 0.05 * diff(zl)], 'Box', 'off', 'FontSize', Q.fontsize - 1, ...
    'TickDir', 'out', 'LineWidth', 0.4, 'XGrid', 'on', 'YGrid', 'on', 'ZGrid', 'on', 'GridAlpha', 0.15);
view(ax, Q.topo_view);
if compact
    % an inset is a few cm wide: the box ends and the scale's ends and zero
    % are all the ticks it has room for
    set(ax, 'XTick', [ceil(min(T.u)) floor(max(T.u))], 'YTick', [ceil(min(T.v)) floor(max(T.v))], ...
        'ZTick', unique([round(clim(1), 1) 0 round(clim(2), 1)]));
    xlabel(ax, [T.axes{1} ' (mm)'], 'FontSize', Q.fontsize - 1);
    ylabel(ax, [T.axes{2} ' (mm)'], 'FontSize', Q.fontsize - 1);
else
    xlabel(ax, axis_label(T.axes{1}), 'FontSize', Q.fontsize - 1);
    ylabel(ax, axis_label(T.axes{2}), 'FontSize', Q.fontsize - 1);
    zlabel(ax, Q.zlabel, 'FontSize', Q.fontsize - 1, 'Interpreter', 'tex');
end

end


function draw_topography_floor(ax, T, zfloor, cm, clim)
% Contours of the surface and its gradient field on the floor.

Z = T.Z;
if nnz(isfinite(Z)) < 9, return, end
levels = linspace(clim(1), clim(2), 13);
C = contourc(T.u, T.v, Z', levels);
k = 1;
while k < size(C, 2)
    lev = C(1, k); n = C(2, k);
    seg = C(:, k + 1:k + n);
    ci = 1 + round((lev - clim(1)) / diff(clim) * (size(cm, 1) - 1));
    ci = min(max(ci, 1), size(cm, 1));
    plot3(ax, seg(1, :), seg(2, :), repmat(zfloor, 1, n), '-', 'Color', cm(ci, :), 'LineWidth', 0.6);
    k = k + n + 1;
end

[gv, gu] = gradient(Z', mean(diff(T.v)), mean(diff(T.u)));     % rows of Z' run along v
step = max(1, round(numel(T.u) / 14));
[UU, VV] = meshgrid(T.u, T.v);
sel = false(size(UU)); sel(1:step:end, 1:step:end) = true;
sel = sel & isfinite(gu) & isfinite(gv);
if any(sel(:))
    quiver3(ax, UU(sel), VV(sel), repmat(zfloor, nnz(sel), 1), gu(sel), gv(sel), zeros(nnz(sel), 1), ...
        0.8, 'Color', [0.1 0.1 0.1], 'LineWidth', 0.4, 'MaxHeadSize', 0.6);
end

end


function draw_topography_inset(fig, pos, T, ~, cm, clim, Q)
% The ROI in 3-D (Kragel et al. 2023 Fig 1c inset), each vertex coloured by
% the FITTED topography at its own projected coordinates - the same values,
% on the same scale, as the surface beside it. Raw voxel values would
% saturate against limits set by the smoothed fit. (render_on_surface is not
% used: it samples the volume at the vertices, which on an ROI's boundary lie
% just outside it, where a masked map is empty.)

ax = axes('Parent', fig, 'Position', pos);
axes(ax);
h = isosurface(to_region(T.roi));
h = h(isgraphics(h));
% only the hemispheres the surface was fitted on. One patch can hold both
% (isosurface joins the two striata), so faces are kept by vertex side.
if any(strcmp(T.hemisphere, {'left', 'right'}))
    sg = 1 - 2 * strcmp(T.hemisphere, 'left');
    for k = 1:numel(h)
        V = get(h(k), 'Vertices'); F = get(h(k), 'Faces');
        set(h(k), 'Faces', F(all(sign(reshape(V(F, 1), size(F))) == sg, 2), :));
    end
end
for k = 1:numel(h)
    V = get(h(k), 'Vertices');
    z = fnval(T.spline, [axis_values(V, T.axes{1})'; axis_values(V, T.axes{2})'])';
    w = (z - clim(1)) / diff(clim);
    ci = 1 + round(min(max(w, 0), 1) * (size(cm, 1) - 1));
    set(h(k), 'FaceVertexCData', cm(ci, :), 'FaceColor', 'interp', 'EdgeColor', 'none', 'FaceAlpha', 1);
end
axis(ax, 'image', 'vis3d', 'off');
view(ax, [0 75]);
lightRestoreSingle(ax);
material(ax, 'dull');
% L / R above each volume, on the figure overlay rather than in this 3-D
% axes: the volumes are in one patch, so the side is taken per vertex
% (vertices still in a face: a one-hemisphere inset keeps the other's vertices)
V = cell2mat(arrayfun(@(x) used_vertices(x), h(:), 'UniformOutput', false));
overlay = findall(fig, 'Type', 'axes', 'Tag', 'LaBGAScore_overlay');
if isempty(V) || isempty(overlay), return, end
drawnow;
px = matlab.graphics.chart.internal.convertDataSpaceCoordsToViewerCoords(ax, V');
u = get(fig, 'Units'); set(fig, 'Units', 'pixels'); fp = get(fig, 'Position'); set(fig, 'Units', u);
pt = [px(1, :)' / fp(3), px(2, :)' / fp(4)];
for sg = [-1 1]
    k = sign(V(:, 1)) == sg;
    if ~any(k), continue, end
    text(overlay(1), mean(pt(k, 1)), max(pt(k, 2)) + 0.01, char('L' * (sg < 0) + 'R' * (sg > 0)), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', Q.fontsize - 1);
end

end



function V = used_vertices(h)
% isosurface returns numeric handles here, so no dot indexing

V = get(h, 'Vertices');
F = get(h, 'Faces');
V = V(unique(F(:)), :);

end



%% =======================================================================
%  PRESETS
%  =======================================================================

function [fig, o, info] = draw_preset(obj, P, info)

before = findall(groot, 'Type', 'figure');

if strcmp(P.layout, 'multirow')
    if ~iscell(obj)
        error('LaBGAScore_brain_figure:multirow', 'layout ''multirow'' needs a cell array of maps, one per row.');
    end
    n = numel(obj);
    o = canlab_results_fmridisplay([], 'multirow', n, 'overlay', P.underlay, 'noverbose');
    C = resolve_colour(obj, P);
    for k = 1:n
        layers = build_layers(obj{k}, C, P);
        for i = 1:numel(layers)
            src = layers{i}.src;
            if isempty_layer(src), continue, end
            o = addblobs(o, to_region(src), layers{i}.args{:}, 'wh_montages', [2 * k - 1, 2 * k], 'noverbose');
        end
        if numel(P.titles) >= k
            o = title_montage(o, 2 * k, P.titles{k});
        end
    end
else
    o = canlab_results_fmridisplay([], P.layout, 'overlay', P.underlay, 'noverbose');
    C = resolve_colour(obj, P);
    layers = build_layers(obj, C, P);
    for i = 1:numel(layers)
        src = layers{i}.src;
        if isempty_layer(src), continue, end
        o = addblobs(o, to_region(src), layers{i}.args{:}, 'noverbose');
    end
end

% canlab_results_fmridisplay builds its own fmridisplay, which opens the
% interactive controller as a uifigure (see 'nocontroller' in the main
% function); close it, and keep the montage figure
created = setdiff(findall(groot, 'Type', 'figure'), before);
ctrl = created(arrayfun(@(f) matlab.ui.internal.isUIFigure(f) || ...
    startsWith(get(f, 'Name'), 'CANlab display controller'), created));
fig = setdiff(created, ctrl);
delete(ctrl);
if isempty(fig), fig = gcf; else, fig = fig(end); end
set(fig, 'Visible', onoff(P.visible), 'Color', 'w', 'InvertHardcopy', 'off');

% Resize to the width requested. CanlabCore placed its axes in normalized
% units on a tall canvas and filled a thin strip of it, so the height is then
% CROPPED to the axes that hold content, with room for the colour key below:
% otherwise a compact2 row comes out as a strip in a mostly empty page.
W = figure_width_mm(P.width);
set(fig, 'Units', 'centimeters');
pos = get(fig, 'Position');
H = W * pos(4) / pos(3);
set(fig, 'Position', [2 2 W / 10 H / 10]);

axs = findobj(fig, 'Type', 'axes');
axs = axs(arrayfun(@(a) ~isempty(get(a, 'Children')), axs));
set(axs, 'Units', 'normalized');
ap = cell2mat(get(axs, 'Position'));
% CanlabCore's slice axes are taller than the figure (1.16 of its height for
% compact2) with the slice letterboxed inside, so their Position says nothing
% about where the pixels are. Project each axes' data box instead.
drawnow;
ybox = zeros(numel(axs), 2);
for k = 1:numel(axs)
    a = axs(k);
    [cx, cy, cz] = ndgrid(get(a, 'XLim'), get(a, 'YLim'), get(a, 'ZLim'));
    pts = zeros(8, 2);
    for c = 1:8, pts(c, :) = data_to_fig(a, [cx(c) cy(c) cz(c)], fig); end
    ybox(k, :) = [min(pts(:, 2)) max(pts(:, 2))];
end
ymin = max(min(ybox(:, 1)), 0); ymax = min(max(ybox(:, 2)), 1);
key_h = 0;
if ~strcmpi(P.colorbar, 'none'), key_h = 4 + P.fontsize * 0.9; end
Hn = P.height;
if isempty(Hn), Hn = (ymax - ymin) * H + key_h + 2; end
scale = ((Hn - key_h - 2) / ((ymax - ymin) * H));    % 1 unless 'height' was given
for k = 1:numel(axs)
    q = ap(k, :);
    q(2) = (key_h + 1 + (q(2) - ymin) * H * scale) / Hn;
    q(4) = q(4) * H * scale / Hn;
    set(axs(k), 'Position', q);
end
set(fig, 'Position', [2 2 W / 10 Hn / 10]);
set(findall(fig, 'Type', 'text'), 'FontName', P.fontname, 'FontSize', P.fontsize);

info.colormode = C.mode;
info.cmaprange = C.range;
if isfield(C, 'bayes'), info.bayes = C.bayes; end

% colorbar below, from the same mapping
overlay = make_overlay_axes(fig);
key_rect = [0.3, 1 / Hn, 0.4, max(key_h - 1, 1) / Hn];
if ~strcmpi(P.colorbar, 'none')
    draw_key(fig, overlay, C, layer_legend(struct('layers', {build_layers(first_map(obj), C, P)})), key_rect, 'south', P);
end

end


function m = first_map(obj)

if iscell(obj), m = obj{1}; else, m = obj; end

end



%% =======================================================================
%  EXPORT
%  =======================================================================

function raster_container(fig, fname, fmt, res, rflag, size_cm)
% A pdf / svg page of exactly the figure's size holding the figure as ONE
% image at the requested dpi. print's '-image' option would do this directly
% but does not exist in R2021a; this way also works headless, where there is
% no OpenGL to rasterise with.

tmp = [tempname '.png'];
print(fig, tmp, '-dpng', res, rflag);
img = imread(tmp);
delete(tmp);
f2 = figure('Visible', 'off', 'Color', 'w', 'Units', 'centimeters', 'Position', [2 2 size_cm], ...
    'PaperUnits', 'centimeters', 'PaperPositionMode', 'manual', 'PaperPosition', [0 0 size_cm], ...
    'PaperSize', size_cm, 'InvertHardcopy', 'off');
ax = axes('Parent', f2, 'Position', [0 0 1 1]);
image(ax, img);
set(ax, 'XLim', [0.5 size(img, 2) + 0.5], 'YLim', [0.5 size(img, 1) + 0.5], 'Visible', 'off');
print(f2, fname, ['-d' fmt], '-painters');
close(f2);

end


function svg_size_in_mm(fname, size_mm)
% print writes an svg's width and height in screen PIXELS at the session's
% DPI (958 x 575 for 183 x 110 mm at 133 dpi), which a viewer assuming 96 dpi
% shows at 253 mm. Restate them in mm, keeping the pixel coordinates as the
% viewBox, so the file opens at the size the journal asked for.

txt = fileread(fname);
w = regexp(txt, '<svg[^>]*?\swidth="([\d.]+)(px)?"', 'tokens', 'once');
h = regexp(txt, '<svg[^>]*?\sheight="([\d.]+)(px)?"', 'tokens', 'once');
if isempty(w) || isempty(h), return, end
if isempty(regexp(txt, '<svg[^>]*\sviewBox=', 'once'))
    txt = regexprep(txt, '<svg\s', sprintf('<svg viewBox="0 0 %s %s" ', w{1}, h{1}), 'once');
end
txt = regexprep(txt, '(<svg[^>]*?\s)width="[\d.]+(px)?"', sprintf('$1width="%.2fmm"', size_mm(1)), 'once');
txt = regexprep(txt, '(<svg[^>]*?\s)height="[\d.]+(px)?"', sprintf('$1height="%.2fmm"', size_mm(2)), 'once');
fid = fopen(fname, 'w');
fwrite(fid, txt);
fclose(fid);

end


function lbf_warn(id, varargin)
% A warning that is also recorded in info.warnings, so it survives into the
% record of what was drawn rather than scrolling past in a headless log.

msg = sprintf(varargin{:});
warning(id, '%s', msg);
w = getappdata(groot, 'LaBGAScore_brain_figure_warnings');
w{end+1} = msg;
setappdata(groot, 'LaBGAScore_brain_figure_warnings', w);

end


function info = finish_figure(fig, P, info)

info.warnings = getappdata(groot, 'LaBGAScore_brain_figure_warnings');

% stray empty axes (e.g. created by a gca call inside a drawing routine)
% would show as a box with ticks
for a = findobj(fig, 'Type', 'axes')'
    if isempty(get(a, 'Children')) && strcmp(get(a, 'Visible'), 'on') && ~strcmp(get(a, 'Tag'), 'LaBGAScore_overlay')
        set(a, 'Visible', 'off');
    end
end

set(fig, 'Units', 'centimeters');
pos = get(fig, 'Position');
info.size_mm = pos(3:4) * 10;
info.dpi = P.dpi;
info.renderer = get(fig, 'Renderer');
info.headless = ~feature('ShowFigureWindows');

if isempty(P.save), return, end

base = char(P.save);
[d, ~, ext] = fileparts(base);
if ~isempty(ext), base = base(1:end - numel(ext)); end
if ~isempty(d) && ~exist(d, 'dir'), mkdir(d); end

% exact physical size: print honours PaperPosition, exportgraphics crops to
% the content and so would change the width a journal asked for
set(fig, 'PaperUnits', 'centimeters', 'PaperPositionMode', 'manual', ...
    'PaperPosition', [0 0 pos(3:4)], 'PaperSize', pos(3:4));

% raster through OpenGL where there is one: print otherwise falls back to
% painters, ~2x slower on CANlab slices (which are surface objects). Headless
% MATLAB has no OpenGL at all, and painters is then the only choice.
rflag = '-painters';
try
    gl = opengl('data');
    if ~strcmpi(gl.Renderer, 'None'), rflag = '-opengl'; end
catch
end
info.export_renderer = rflag(2:end);
res = sprintf('-r%d', P.dpi);

for f = 1:numel(P.formats)
    fmt = lower(P.formats{f});
    fname = [base '.' fmt];
    t0 = tic;
    switch fmt
        case 'png', print(fig, fname, '-dpng', res, rflag);
        case 'tif', print(fig, fname, '-dtiff', res, rflag);
        case {'pdf', 'svg'}
            if P.vector
                print(fig, fname, ['-d' fmt], '-painters');
            else
                raster_container(fig, fname, fmt, res, rflag, pos(3:4));
            end
            if strcmp(fmt, 'svg'), svg_size_in_mm(fname, pos(3:4) * 10); end
        otherwise
            error('LaBGAScore_brain_figure:format', 'Unsupported format ''%s''.', fmt);
    end
    info.files{end+1} = fname;
    if P.verbose, fprintf('LaBGAScore_brain_figure: wrote %s (%.0f s)\n', fname, toc(t0)); end
end

end
