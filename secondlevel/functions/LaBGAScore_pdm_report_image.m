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
% Delegated to LaBGAScore_blob_montage, which the SVM scripts also use, so the
% montage styling, font scaling, figure tagging and sizing live in ONE place.
% 'noregioncenters' here: the regioncenters montage comes after the table below.
LaBGAScore_blob_montage(d, region(d), sprintf('%s (%d voxels)', label, nsig), ...
    'fontscale', fontscale, 'max_regioncenters', max_regioncenters, 'noregioncenters');

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
% Same delegation; 'regioncentersonly' skips the overview montage already drawn
% above, so the reporting order stays montage -> table -> regioncenters.
LaBGAScore_blob_montage(d, r, label, 'fontscale', fontscale, ...
    'max_regioncenters', max_regioncenters, 'regioncentersonly');

end
