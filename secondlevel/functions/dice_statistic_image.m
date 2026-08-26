function [dice, info] = dice_statistic_image(img1, img2, varargin)
% dice_statistic_image  Dice overlap between two thresholded statistic images.
%
%   [dice, info] = dice_statistic_image(img1, img2)
%   [dice, info] = dice_statistic_image(img1, img2, 'p_threshold', 0.01)
%   [dice, info] = dice_statistic_image(img1, img2, 'binary_mask', {mA, mB})
%
%   Computes the Dice similarity coefficient between two CANlab
%   statistic_image objects, either by thresholding their p-values or from
%   explicitly supplied binary masks.
%
%       dice = 2 * |A ∩ B| / (|A| + |B|)
%
%   INPUTS
%     img1, img2  statistic_image objects sharing one voxel space
%
%   OPTIONAL NAME-VALUE PAIRS
%     'p_threshold'  p-value threshold, default 0.05. Voxels with p BELOW the
%                    threshold count as suprathreshold; the comparison is
%                    strict (p < threshold), matching
%                    thresholded_fmri_data_from_statistic_image and
%                    dice_statistic_image_by_roi.
%     'binary_mask'  {maskA, maskB} logical vectors, which override p-value
%                    thresholding entirely.
%
%   OUTPUTS
%     dice  Dice coefficient, NaN when neither image has a suprathreshold
%           voxel (the coefficient is undefined there rather than zero)
%     info  struct with .nA, .nB, .nOverlap, .maskA, .maskB, .p_threshold
%
%   NaN p-values are treated as non-significant, since any comparison against
%   NaN is false.
%
%   See also dice_statistic_image_by_roi, thresholded_fmri_data_from_statistic_image.

%% ---------------- defaults ----------------
p_thresh = 0.05;
binary_mask = [];

for i = 1:2:numel(varargin)
   switch lower(varargin{i})
       case 'p_threshold'
           p_thresh = varargin{i+1};
       case 'binary_mask'
           binary_mask = varargin{i+1};
       otherwise
           error('dice_statistic_image:unknownOption', ...
               'Unknown option: %s', varargin{i});
   end
end

%% ---------------- checks ----------------
if ~isa(img1, 'statistic_image') || ~isa(img2, 'statistic_image')
   error('dice_statistic_image:inputType', ...
       'Inputs must be CANlab statistic_image objects.');
end

if numel(img1.dat) ~= numel(img2.dat)
   error('dice_statistic_image:imageSizeMismatch', ...
       'img1 has %d voxels but img2 has %d; they must share a voxel space.', ...
       numel(img1.dat), numel(img2.dat));
end

%% ---------------- build binary masks ----------------
if ~isempty(binary_mask)
   % User-supplied masks
   if ~iscell(binary_mask) || numel(binary_mask) ~= 2
       error('dice_statistic_image:binaryMask', ...
           'binary_mask must be a cell array with two elements.');
   end

   maskA = logical(binary_mask{1}(:));
   maskB = logical(binary_mask{2}(:));

   if numel(maskA) ~= numel(img1.dat) || numel(maskB) ~= numel(img2.dat)
       error('dice_statistic_image:binaryMaskSize', ...
           ['Supplied masks have %d and %d elements but the images have %d ' ...
            'voxels each.'], numel(maskA), numel(maskB), numel(img1.dat));
   end

else
   % Threshold using p-values. NaN < threshold is false, so non-finite
   % p-values fall out as non-significant without needing a separate mask.
   if isempty(img1.p) || isempty(img2.p)
       error('dice_statistic_image:noPvalues', ...
           'Images carry no p-values; supply ''binary_mask'' explicitly.');
   end

   maskA = img1.p(:) < p_thresh;
   maskB = img2.p(:) < p_thresh;
end

%% ---------------- Dice coefficient ----------------
nA = sum(maskA);
nB = sum(maskB);
nOverlap = sum(maskA & maskB);

if nA + nB == 0
   dice = NaN; % undefined if neither map has suprathreshold voxels
else
   dice = 2 * nOverlap / (nA + nB);
end

%% ---------------- output info ----------------
info = struct();
info.nA = nA;
info.nB = nB;
info.nOverlap = nOverlap;
info.maskA = maskA;
info.maskB = maskB;
info.p_threshold = p_thresh;

end
