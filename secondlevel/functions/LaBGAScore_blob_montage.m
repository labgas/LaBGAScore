function LaBGAScore_blob_montage(dat, r, label, varargin)
% LaBGAScore_blob_montage  Overview + regioncenters montages for a thresholded map.
%
% Display only: it draws, it does not threshold, label or tabulate. Callers that
% already built a region object and table (thresholded_fmri_data_from_pval_nii,
% LaBGAScore_pdm_report_image) pass them in, so nothing is recomputed and no
% table is printed twice.
%
% Nothing is drawn when there are no suprathreshold regions, so adding a call to
% an analysis that comes out null costs nothing and leaves no empty figures.
%
% :Usage:
% ::
%     LaBGAScore_blob_montage(fmri_dat, region_obj, 'AUC p_unc < 0.001')
%
% :Inputs:
%   **dat:**    thresholded fmri_data (or fmri_data_st) object
%   **r:**      region object built from dat; [] or empty skips everything
%   **label:**  string used in the montage titles
%
% :Optional inputs:
%   **'max_regioncenters'**  default 21. Above this the regioncenters montage is
%                            skipped: one titled panel per region stops being
%                            readable, and stops being quick.
%   **'fontscale'**          default 2/3. Montage titles are tuned for a larger
%                            canvas than headless publishing provides.
%   **'noregioncenters'**    overview montage only.
%
% -------------------------------------------------------------------------
% by: Lukas Van Oudenhove  |  KU Leuven, September 2026
% LaBGAScore_blob_montage.m   v1.0   last modified: 2026/09/17

max_regioncenters = 21; fontscale = 2/3; doregioncenters = true;
for i = 1:numel(varargin)
    if ~ischar(varargin{i}) && ~isstring(varargin{i}), continue, end
    switch lower(char(varargin{i}))
        case 'max_regioncenters', max_regioncenters = varargin{i+1};
        case 'fontscale',         fontscale = varargin{i+1};
        case 'noregioncenters',   doregioncenters = false;
    end
end

if isempty(dat) || isempty(r), return, end
if isa(dat,'fmri_data') || isa(dat,'fmri_data_st')
    if sum(dat.dat ~= 0) == 0, return, end
end

splitc = {[.1 .8 .8] [.1 .1 .8] [.9 .4 0] [1 1 0]};

% ---- overview montage ---------------------------------------------------
o2 = canlab_results_fmridisplay([], 'compact');
o2 = addblobs(o2, r, 'splitcolor', splitc);
[o2, th] = title_montage(o2, 5, label); %#ok<ASGLU>
if ~isempty(th) && all(ishandle(th))
    set(th, 'FontSize', get(th(1),'FontSize') * fontscale);
end
set(gcf, 'Tag', [matlab.lang.makeValidName(label) '_montage']);
plugin_set_figure_size;
drawnow, snapnow

% ---- regioncenters montage ---------------------------------------------
if ~doregioncenters, return, end

if numel(r) < max_regioncenters
    fprintf('\n  MONTAGE REGIONCENTERS: %s, %d regions\n\n', label, numel(r));
    o3 = montage(r, 'regioncenters', 'splitcolor', splitc); %#ok<NASGU>
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
