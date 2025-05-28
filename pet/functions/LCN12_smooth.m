function smooth_image = LCN12_smooth(image,smooth_kernel,mask)
% LCN12_smooth.m
% 
% This function will smooth an image with a 3D Gaussian kernel optionally 
% within a mask
% INPUT:
%   image = 3D matrix
%   smooth_kernel = array of Gaussian kernel size (in voxels)
%   Optional input: mask = 3D binary mask
% OUTPUT:
%   smooth_image = smoothed image of the same size as image
%
% Important:
%  - requires the installation of SPM12 - http://www.fil.ion.ucl.ac.uk/spm/
%  - the path of SPM and this routine should be included in the Matlab path
%
% The package contains software (SPM12) developed under the auspices of The
% Wellcome Department of Imaging Neuroscience, a department of the
% Institute of Neurology at University College London. The copyright of
% this software remains with that of SPM12, see
% http://www.fil.ion.ucl.ac.uk/spm/.
%    
% This routine is supplied as is. 
% Comments or questions can be send to:
% Patrick.Dupont@med.kuleuven.be
%
%__________________________________________________________________________
%
% author: 	Patrick Dupont
%           Labo for cognitive neurology, KU LEUVEN 
% date: 	April, 2021
% history: 	June 2021: bug fix if no mask was given as argument
%__________________________________________________________________________
% @(#)LCN12_smooth.m	   v0.11                  last modified: 2021/06/16

if nargin == 3
   % smooth within the mask
   % this can be done using the smooth of the original image (assuming
   % zeros outside the mask) and then divide this by the value obtained by
   % smoothing the mask itself (assuming binary mask)
   % This idea is coming from Mark Jenkinson in the FSL mailing list
   mask(mask>0) = 1;
   mask(mask <0) = 0;
   mask = double(mask);
   imagenew = image.*(mask == 1);
   smooth_imagenew = zeros(size(imagenew));
   smooth_mask     = zeros(size(mask));
   spm_smooth(imagenew,smooth_imagenew,smooth_kernel);   
   spm_smooth(mask,smooth_mask,smooth_kernel);   
   smooth_image = zeros(size(image));
   smooth_image(mask == 1) = smooth_imagenew(mask == 1)./smooth_mask(mask == 1); 
else
   smooth_image = zeros(size(image));
   spm_smooth(image,smooth_image,smooth_kernel);   
end
end