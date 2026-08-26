function roi_table = dice_statistic_image_by_roi(imgA, imgB, roi_img, varargin)
% dice_statistic_image_by_roi  Per-ROI Dice overlap between two statistic images.
%
% Compute per-ROI Dice coefficients between two statistic_image objects.
%
% INPUTS
%  imgA, imgB : statistic_image
%  roi_img   : ROI definition (statistic_image or fmri_data)
%
% OPTIONAL NAME–VALUE
%  'binary_mask' : {maskA, maskB}
%  'p_threshold' : p-value threshold (default = 0.05)
%
% OUTPUT
%  roi_table : table with columns:
%      ROI      ROI_label
%      Dice     Dice coefficient
%      nA       Voxels in A within ROI
%      nB       Voxels in B within ROI
%      nOverlap Overlapping voxels in ROI
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
   maskA = imgA.p < p_thresh;
   maskB = imgB.p < p_thresh;
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
roi_labels = unique(roi_img.dat);
roi_labels(roi_labels == 0 | isnan(roi_labels)) = [];

nROI     = numel(roi_labels);
ROI      = zeros(nROI,1);
Dice     = nan(nROI,1);
nA       = zeros(nROI,1);
nB       = zeros(nROI,1);
nOverlap = zeros(nROI,1);

for r = 1:numel(roi_labels)
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

