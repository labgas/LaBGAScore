function p_fwe = tfce_fwe_from_null(tfce_real, TFCE_null_max)
% tfce_fwe_from_null  Voxelwise FWE-corrected p from a max-statistic null.
%
%   p_fwe = tfce_fwe_from_null(tfce_real, TFCE_null_max)
%
%   Implements the standard permutation FWE correction (Smith & Nichols 2009):
%   each voxel's observed statistic is compared with the distribution of the
%   PER-PERMUTATION MAXIMUM statistic, so a voxel is significant only if it
%   beats what the strongest voxel anywhere in the brain achieved by chance.
%   That is what makes the correction familywise over the whole search volume.
%
%       p_fwe(v) = ( #{ TFCE_null_max >= tfce_real(v) } + 1 ) / (nPerm + 1)
%
%   The +1 in both terms is the standard Phipson & Smyth (2010) correction: the
%   observed data are themselves one realization under the null, so p is never
%   exactly zero and the test stays valid at the smallest attainable p, 1/(nPerm+1).
%
%   INPUTS
%     tfce_real      nVox x 1, observed TFCE statistic
%     TFCE_null_max  nPerm x 1, maximum TFCE statistic within each permutation
%
%   OUTPUT
%     p_fwe          nVox x 1, FWE-corrected p-values, in the voxel order of
%                    tfce_real. ALREADY CORRECTED - threshold these directly
%                    (e.g. at .05); running FDR on top would correct twice.
%
%   The p resolution is 1/(nPerm+1), so nPerm must be at least a few hundred for
%   a .05 threshold to be meaningful, and 5000-10000 is usual for final inference.
%
%   Computed by sorting rather than by comparing every voxel against every
%   permutation: an nPerm x nVox logical array is 9.4 GB at 10,000 permutations
%   and 235,807 voxels, whereas this is O(n log n) in a few hundred MB.
%
%   Factored out of group_tfce_from_subject_maps so that the same arithmetic can
%   be applied retrospectively to TFCE results saved before FWE was implemented -
%   those store TFCE_real and TFCE_null_max, which is all this needs.
%
%   See also group_tfce_from_subject_maps.

% AUTHOR
% Lukas Van Oudenhove, KU Leuven
%
% ..

nPerm = numel(TFCE_null_max);

[real_sorted, real_order] = sort(double(tfce_real(:)));
null_sorted = sort(double(TFCE_null_max(:)));

% histcounts over edges [-inf; sorted real values; inf]: the cumulative count up
% to bin i is the number of null maxima strictly below the i-th real value
bin_counts  = histcounts(null_sorted, [-inf; real_sorted; inf]);
n_null_below = cumsum(bin_counts(1:end-1))';

p_fwe = zeros(size(real_sorted));
p_fwe(real_order) = (nPerm - n_null_below + 1) ./ (nPerm + 1);

end % function
