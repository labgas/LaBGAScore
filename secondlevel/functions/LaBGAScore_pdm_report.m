function out = LaBGAScore_pdm_report(pdm, dat_template, condir, varargin)
% LaBGAScore_pdm_report  Report thresholded PDM maps the way c2a reports GLM maps.
%
% The PDM toolbox writes weight images and little else: plotPDM gives a bare
% montage with no labelling, no region table, and no per-region detail. Second
% level GLM results are reported through canlab_results_fmridisplay, a labelled
% region table, and a regioncenters montage. This brings PDM output up to the
% same standard so the two can be read side by side.
%
% :Usage:
% ::
%     out = LaBGAScore_pdm_report(pdm, dat_template, condir, 'atlas', atlas_obj, ...)
%
% :Inputs:
%   **pdm:**           the struct returned by multivariateMediation, carrying
%                      Wfull, boot.p and pThreshold
%   **dat_template:**  an fmri_data object in the analysis space, used to turn
%                      the weight vectors back into images
%   **condir:**        directory to write PDM<k>.nii into
%
% :Optional inputs:
%   **'atlas'**             atlas object used to label regions. Default: none,
%                           in which case the region table falls back to its
%                           own default atlas.
%   **'k'**                 cluster extent threshold in voxels, default 0
%   **'max_regioncenters'** regioncenters montages stop being readable, and
%                           stop being quick, once there are many regions.
%                           Default 21, matching c2a.
%   **'titlestr'**          string identifying the analysis in figure titles
%   **'fontscale'**         multiplier on montage title font size, default 2/3,
%                           because the default is tuned for a larger canvas
%                           than headless publishing provides
%   **'nowrite'**           do not write the .nii files (they already exist)
%
% :Outputs:
%   **out:**  struct array, one element per PDM, with fields region, table,
%             table_cov, nsig and fname
%
% -------------------------------------------------------------------------
% by: Lukas Van Oudenhove
% date: KU Leuven, September 2026
% -------------------------------------------------------------------------
% LaBGAScore_pdm_report.m         v1.0
% last modified: 2026/09/17

atlas_obj = [];
k_threshold = 0;
max_regioncenters = 21;
titlestr = '';
fontscale = 2/3;
dowrite = true;

for i = 1:2:numel(varargin)
    switch lower(varargin{i})
        case 'atlas',             atlas_obj = varargin{i+1};
        case 'k',                 k_threshold = varargin{i+1};
        case 'max_regioncenters', max_regioncenters = varargin{i+1};
        case 'titlestr',          titlestr = varargin{i+1};
        case 'fontscale',         fontscale = varargin{i+1};
        case 'nowrite',           dowrite = false; i = i - 1; %#ok<FXSET>
        otherwise
            error('unrecognised option ''%s''', varargin{i});
    end
end

out = struct('region',{},'table',{},'table_cov',{},'nsig',{},'fname',{});

if ~isfield(pdm,'boot') || ~isfield(pdm.boot,'p')
    fprintf('\n  no bootstrap performed, so no thresholded PDM images or tables\n');
    return
end

for k = 1:numel(pdm.boot.p)

    d = dat_template;
    d.dat = pdm.Wfull{k} .* (pdm.boot.p{k} < pdm.pThreshold(k));
    nsig = sum(d.dat ~= 0);

    fprintf('\n\n');
    fprintf('==== PDM%d: %d voxel(s) below the bootstrap threshold p < %.4g ====\n', ...
        k, nsig, pdm.pThreshold(k));
    fprintf('\n');

    fname = fullfile(condir, sprintf('PDM%d.nii', k));
    if dowrite
        write(d, 'fname', fname, 'overwrite');
    end

    out(k).nsig = nsig; %#ok<AGROW>
    out(k).fname = fname; %#ok<AGROW>
    out(k).region = []; out(k).table = []; out(k).table_cov = [];

    if nsig == 0
        fprintf('  nothing survives threshold, so no montage or table for this PDM\n');
        continue
    end

    thistitle = sprintf('PDM%d %s (p < %.4g)', k, titlestr, pdm.pThreshold(k));

    % Reporting itself lives in LaBGAScore_pdm_report_image so the identical
    % montage/table/regioncenters output can be regenerated later from the
    % written PDM*.nii without refitting the bootstrap.
    rep = LaBGAScore_pdm_report_image(d, thistitle, 'atlas', atlas_obj, ...
        'k', k_threshold, 'max_regioncenters', max_regioncenters, 'fontscale', fontscale);

    out(k).region = rep.region; out(k).table = rep.table; out(k).table_cov = rep.table_cov;

end

end

