function [img,V] = LCN12_read_image(filename,Vref)
% LCN12_read_image
%
% this script will read an image in matlab.
%
% INPUT
%           filename = name of the file (either full path or name if you
%                      are in the correct directory). If the file is a 4D
%                      file, you have to specify the volume you want to
%                      read in the filename (see SPM12 documentation).
%           Vref = reference structure obtained using spm_vol and if you
%                  want that the image is resampled according to the 
%                  reference image. This is an optional argument
%
% OUTPUT
%           img = 3D matrix of an image
%           V   = structure from SPM12 with the information of the mapped
%                 volume. If you have specified a reference structure Vref, 
%                 this will return an empty value.
%
% Important:
%  - requires the installation of SPM12 - http://www.fil.ion.ucl.ac.uk/spm/
%  - the path of SPM and this routine should be included in the Matlab path
%
% The package contains software (SPM12) developed under the auspices of The
% Wellcome Department of Imaging Neuroscience, a department of the
% Institute of Neurology at University College London. The copyright of
% this software remains with that of SPM, see http://www.fil.ion.ucl.ac.uk/spm/.
%    
% This routine is supplied as is. 
%
% Comments or questions can be send to:
% Patrick.Dupont@med.kuleuven.be
%__________________________________________________________________________
%
% author: Patrick Dupont
% date:   October 17, 2015
% history: 
%__________________________________________________________________________
% @(#)LCN12_read_image.m           v0.0           last modified: 2015/10/17

V = spm_vol(filename);
if nargin == 2
   nr_planes = Vref.dim(3);
   % initialize
   img = zeros(Vref.dim(1),Vref.dim(2),nr_planes);
   for z = 1:nr_planes
       B = spm_matrix([0 0 -z 0 0 0 1 1 1]);
       M = inv(B*inv(Vref.mat)*V.mat);
       d = spm_slice_vol(V,M,Vref.dim(1:2),1);
       img1 = reshape(d,Vref.dim(1:2));
       img(:,:,z) = img1;
   end
   V = [];
else
   img = spm_read_vols(V);
end