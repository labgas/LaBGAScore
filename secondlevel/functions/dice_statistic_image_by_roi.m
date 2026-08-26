function roi_table = dice_statistic_image_by_roi(imgA, imgB, roi_img, varargin)
% dice_statistic_image_by_roi  Per-ROI Dice overlap between two statistic images.
%
%   roi_table = dice_statistic_image_by_roi(imgA, imgB, roi_img)
%   roi_table = dice_statistic_image_by_roi(imgA, imgB, roi_img, 'p_threshold', 0.01)
%   roi_table = dice_statistic_image_by_roi(imgA, imgB, roi_img, 'binary_mask', {mA, mB})
%
%   Computes the Dice similarity coefficient separately within each ROI of a
%   label image. The whole-image counterpart is dice_statistic_image.
%
%   INPUTS
%     imgA, imgB  statistic_image objects sharing one voxel space
%     roi_img     LABEL image in that same space (atlas, statistic_image or
%                 fmri_data). Its .dat must hold integer region labels, with 0
%                 for background; one table row is produced per distinct
%                 non-zero label. Do not pass a continuous statistic map: every
%                 distinct value would be treated as its own ROI.
%
%   OPTIONAL NAME-VALUE PAIRS
%     'p_threshold'  p-value threshold, default 0.05. The comparison is strict
%                    (p < threshold), matching dice_statistic_image and
%                    thresholded_fmri_data_from_statistic_image.
%     'binary_mask'  {maskA, maskB} logical vectors, which override p-value
%                    thresholding entirely.
%
%   OUTPUT
%     roi_table  one row per ROI, with columns
%        ROI       the integer label
%        Dice      Dice coefficient, NaN where neither image has a
%                  suprathreshold voxel in that ROI (undefined, not zero)
%        nA        suprathreshold voxels of imgA within the ROI
%        nB        suprathreshold voxels of imgB within the ROI
%        nOverlap  voxels suprathreshold in both, within the ROI
%
%   All three images are checked to share a voxel space. NaN p-values are
%   treated as non-significant, since any comparison against NaN is false.
%
%   See also dice_statistic_image, thresholded_fmri_data_from_statistic_image.
%

%% defaults
p_thresh = 0.05;
binary_mask = [];

for i = 1:2:numel(varargin)
   switch lower(varargin{i})
       case 'p_threshold'
           p_thresh = varargin{i+1};
       case 'binary_mask'
           binary_mask = varargin{i+1};
       otherwise
           % The sibling dice_statistic_image errors here; this one used to
           % fall through, so a mistyped option name was silently ignored and
           % the default was used instead.
           error('dice_statistic_image_by_roi:unknownOption', ...
               'Unknown option: %s', varargin{i});
   end
end

%% checks
% dice_statistic_image verifies its two inputs share a voxel space; this one
% checked nothing, including that roi_img lives in the same space. A mismatch
% would either error deep inside the ROI loop or, worse, silently compare
% voxels that do not correspond to each other.
if ~isa(imgA,'statistic_image') || ~isa(imgB,'statistic_image')
   error('dice_statistic_image_by_roi:inputType', ...
       'imgA and imgB must be CANlab statistic_image objects.');
end

if numel(imgA.dat) ~= numel(imgB.dat)
   error('dice_statistic_image_by_roi:imageSizeMismatch', ...
       'imgA has %d voxels but imgB has %d; they must share a voxel space.', ...
       numel(imgA.dat), numel(imgB.dat));
end

if numel(roi_img.dat) ~= numel(imgA.dat)
   error('dice_statistic_image_by_roi:roiSizeMismatch', ...
       ['roi_img has %d voxels but the statistic images have %d. Resample the ' ...
        'ROI image into the same space before calling.'], ...
       numel(roi_img.dat), numel(imgA.dat));
end

%% masks
if isempty(binary_mask)
   if isempty(imgA.p) || isempty(imgB.p)
       error('dice_statistic_image_by_roi:noPvalues', ...
           'Images carry no p-values; supply ''binary_mask'' explicitly.');
   end
   % NaN < thresh is false, so non-finite p-values count as non-significant
   maskA = imgA.p(:) < p_thresh;
   maskB = imgB.p(:) < p_thresh;
else
   if numel(binary_mask) ~= 2
       error('dice_statistic_image_by_roi:binaryMask', ...
           'binary_mask must be a cell array with two elements.');
   end
   maskA = logical(binary_mask{1}(:));
   maskB = logical(binary_mask{2}(:));
   if numel(maskA) ~= numel(imgA.dat) || numel(maskB) ~= numel(imgB.dat)
       error('dice_statistic_image_by_roi:binaryMaskSize', ...
           'Supplied masks do not match the image voxel count.');
   end
end

%% ROI labels
roi_labels = unique(roi_img.dat(:));
roi_labels(roi_labels == 0 | isnan(roi_labels)) = [];

% A continuous map passed as roi_img would yield one "ROI" per distinct voxel
% value, silently producing a table with thousands of meaningless rows.
if any(roi_labels ~= round(roi_labels))
    error('dice_statistic_image_by_roi:nonIntegerLabels', ...
        ['roi_img holds non-integer values, so it is not a label image. Pass an ' ...
         'atlas or ROI mask whose .dat contains integer region labels.']);
end

if numel(roi_labels) > numel(roi_img.dat)/10
    warning('dice_statistic_image_by_roi:manyLabels', ...
        ['roi_img has %d distinct labels across %d voxels. If this is not really ' ...
         'a label image, the table below will be meaningless.'], ...
        numel(roi_labels), numel(roi_img.dat));
end

nROI     = numel(roi_labels);
ROI      = zeros(nROI,1);
Dice     = nan(nROI,1);
nA       = zeros(nROI,1);
nB       = zeros(nROI,1);
nOverlap = zeros(nROI,1);

for r = 1:nROI
   lab = roi_labels(r);
   roi_mask = roi_img.dat == lab;

   A = maskA & roi_mask;
   B = maskB & roi_mask;

   nA_r = sum(A);
   nB_r = sum(B);
   nOv = sum(A & B);

   if nA_r + nB_r == 0
       dice_r = NaN;
   else
       dice_r = 2 * nOv / (nA_r + nB_r);
   end

   ROI(r,1) = lab;
   Dice(r,1) = dice_r;
   nA(r,1) = nA_r;
   nB(r,1) = nB_r;
   nOverlap(r,1) = nOv;
end

roi_table = table(ROI, Dice, nA, nB, nOverlap);

end

