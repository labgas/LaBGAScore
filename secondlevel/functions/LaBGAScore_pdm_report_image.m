function out = LaBGAScore_pdm_report_image(d, label, varargin)
% LaBGAScore_pdm_report_image  c2a-style reporting for one thresholded brain map.
%
% Produces the overview montage, the labelled region table, and the
% regioncenters montage that second-level GLM results get from c2a, for a map
% that is ALREADY thresholded. Split out from LaBGAScore_pdm_report so reports
% can be regenerated from written PDM*.nii without refitting the bootstrap,
% which takes over an hour per analysis.
%
% :Usage:
% ::
%     out = LaBGAScore_pdm_report_image(d, 'PDM1 discoverie GM', 'atlas', atl)
%
% :Inputs:
%   **d:**      a thresholded fmri_data object (zeros where not significant)
%   **label:**  string identifying the map, used in figure titles
%
% :Optional inputs:  'atlas', 'k', 'max_regioncenters', 'fontscale'
%                    (see LaBGAScore_pdm_report)
%
% :Output:
%   **out:**  struct with fields region, table, table_cov, nsig
%
% -------------------------------------------------------------------------
% by: Lukas Van Oudenhove   |   KU Leuven, September 2026
% LaBGAScore_pdm_report_image.m   v1.0   last modified: 2026/09/17

atlas_obj = []; k_threshold = 0; max_regioncenters = 21; fontscale = 2/3;
for i = 1:2:numel(varargin)
    switch lower(varargin{i})
        case 'atlas',             atlas_obj = varargin{i+1};
        case 'k',                 k_threshold = varargin{i+1};
        case 'max_regioncenters', max_regioncenters = varargin{i+1};
        case 'fontscale',         fontscale = varargin{i+1};
        otherwise, error('unrecognised option ''%s''', varargin{i});
    end
end

out = struct('region',[],'table',[],'table_cov',[],'nsig',0);

nsig = sum(d.dat ~= 0);
out.nsig = nsig;

fprintf('\n\n');
fprintf('==== %s: %d suprathreshold voxel(s) ====\n', label, nsig);
fprintf('\n');

if nsig == 0
    fprintf('  nothing survives threshold, so no montage or table\n');
    return
end

% ---- overview montage ---------------------------------------------------
% PDM weights are signed. The sign of a whole PDM is arbitrary, but the
% RELATIVE sign across voxels is meaningful, so split the colour map rather
% than showing magnitude only.
o2 = canlab_results_fmridisplay([], 'compact');
o2 = addblobs(o2, region(d), 'splitcolor', {[.1 .8 .8] [.1 .1 .8] [.9 .4 0] [1 1 0]});
[o2, th] = title_montage(o2, 5, sprintf('%s (%d voxels)', label, nsig)); %#ok<ASGLU>
if all(ishandle(th)), set(th, 'FontSize', get(th(1),'FontSize') * fontscale); end
set(gcf, 'Tag', [matlab.lang.makeValidName(label) '_montage']);
plugin_set_figure_size;
drawnow, snapnow

% ---- region table -------------------------------------------------------
r = region(d);
if ~isempty(r) && k_threshold > 0
    r(cat(1, r.numVox) < k_threshold) = [];
end
if isempty(r)
    fprintf('  no regions left after the extent threshold\n');
    return
end

if ~isempty(atlas_obj)
    [rpos, rneg, r_table] = LaBGAScore_region_table_safe(@table, r, 'atlas_obj', atlas_obj);
else
    [rpos, rneg, r_table] = LaBGAScore_region_table_safe(@table, r);
end
r = [rpos rneg];

r_table_cov = [];
if ~isempty(atlas_obj)
    try
        [~, ~, ~, ~, ~, ~, r_table_cov] = table_of_atlas_regions_covered(d, 'atlas', atlas_obj);
    catch ME
        fprintf('  table_of_atlas_regions_covered failed (%s); continuing\n', ME.message);
    end
end

out.region = r; out.table = r_table; out.table_cov = r_table_cov;

% ---- regioncenters montage ---------------------------------------------
% One titled panel per region: unreadable and slow once there are many, so
% gate it exactly as c2a does.
if numel(r) < max_regioncenters
    fprintf('\n  MONTAGE REGIONCENTERS, %s, %d regions\n\n', label, numel(r));
    o3 = montage(r, 'regioncenters', 'splitcolor', {[.1 .8 .8] [.1 .1 .8] [.9 .4 0] [1 1 0]}); %#ok<NASGU>
    % regioncenter panel titles are the worst offenders at the headless canvas
    % size, so scale them too.
    tt = findobj(gcf, 'Type', 'text');
    for z = 1:numel(tt)
        try, set(tt(z), 'FontSize', get(tt(z),'FontSize') * fontscale); catch, end %#ok<CTCH>
    end
    set(gcf, 'Tag', [matlab.lang.makeValidName(label) '_regioncenters']);
    plugin_set_figure_size;
    drawnow, snapnow
else
    fprintf('\n  regioncenters montage skipped: %d regions, at or above the display limit of %d\n\n', ...
        numel(r), max_regioncenters);
end

end
