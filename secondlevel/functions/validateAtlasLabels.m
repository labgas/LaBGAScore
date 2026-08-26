function labels = validateAtlasLabels(atlasData, p, caller)
% validateAtlasLabels  Check that atlas label i really is feature i.
%
%   labels = validateAtlasLabels(atlasData, p, caller)
%
%   The three diagnostic plotters paint results into a NIfTI with
%
%       for i = 1:p
%           mask = atlasData == i;
%           map(mask) = value(i);
%       end
%
%   so the correspondence "atlas integer label i = column i of X" must hold
%   EXACTLY. This function enforces it.
%
%   Inputs
%     atlasData  3D atlas volume as read by spm_read_vols
%     p          number of features (columns of X)
%     caller     char, used in the error/warning identifier and message
%
%   Output
%     labels     the sorted non-zero, non-NaN labels found in the volume
%
%   Behaviour
%     - labels missing from 1:p  -> ERROR. Those features would be painted onto
%       nothing, and every later feature would land on the wrong region.
%     - labels beyond p          -> warning. Those atlas regions simply are not
%       represented in X and stay empty in the map, which is usually intended.
%
%   Why the previous check was not enough
%     It tested only `max(labels) < p`. That misses the case that actually
%     bites: labels with GAPS. With labels [1 2 5 7 9] and p = 5, max(labels)
%     is 9, so no warning fired -- yet features 3 and 4 matched no voxels and
%     were silently dropped, and feature 5 was painted onto atlas label 5,
%     which is the THIRD ROI in the set. The output looked like a normal brain
%     map. Gaps arise whenever an ROI ends up with no voxels after masking or
%     resampling, which is routine.
%
%   See also plot_PLSR_diagnostics_neuroimaging, plot_PLSDA_diagnostics_neuroimaging,
%   plot_ENet_diagnostics_neuroimaging.

if nargin < 3 || isempty(caller)
    caller = 'validateAtlasLabels';
end

labels = unique(atlasData(:));
labels(labels == 0 | isnan(labels)) = [];
labels = sort(labels(:))';

if isempty(labels)
    error([caller ':atlasEmpty'], ...
        'The atlas volume contains no non-zero labels, so results cannot be mapped to regions.');
end

missing = setdiff(1:p, labels);
extra   = setdiff(labels, 1:p);

if ~isempty(missing)
    error([caller ':atlasLabelsNotContiguous'], ...
        ['Atlas labels do not match the feature columns. Labels %s are absent from the ' ...
         'atlas volume, but X has %d features and the mapping assumes atlas label i IS ' ...
         'feature i.\n' ...
         'Features at the missing labels would be dropped, and every feature after the ' ...
         'first gap would be painted onto the wrong region.\n' ...
         'Atlas labels present: %s\n' ...
         'Rebuild the ROI atlas so its labels run 1..%d without gaps (a gap usually means ' ...
         'an ROI lost all of its voxels during masking or resampling), or subset X to match.'], ...
        mat2str(missing), p, mat2str(labels), p);
end

if ~isempty(extra)
    warning([caller ':atlasExtraLabels'], ...
        ['The atlas contains %d label(s) beyond the %d features in X (%s). Those regions ' ...
         'are not represented in the results and will be left empty in the maps.'], ...
        numel(extra), p, mat2str(extra));
end

end
