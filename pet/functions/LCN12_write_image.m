function Vout = LCN12_write_image(img,outputfilename,description,datatype,Vref)
% LCN12_write_image
%
% this script will read an image in matlab.
%
% syntax: Vout = LCN12_write_image(img,outputfilename,description,datatype,Vref)
%
% INPUT
%           img = 3D or 4D matrix of an image
%           outputfilename = name of the outputfile. Full path required 
%                            otherwise it will be written to the current 
%                            working directory. Extension is assumed to be
%                            .nii or .img (but .nii is prefered).
%           description = short string which will be written to the output
%           file.
%           datatype = 2 (uint8), 4 (int16), 8 (int32), 16 (float32), 
%                      64 (float64), .... (see spm_type)
%           Vref = reference structure obtained using spm_vol. This is an 
%                  optional argument
%
% OUTPUT
%           Vout = structure from SPM12 with the information of the mapped
%                  volume. If you have specified a reference structure Vref, 
%                  this will return an empty value.
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
% history: Feb 2019: also possible to write 4D nifti files
%__________________________________________________________________________
% @(#)LCN12_write_image.m           v0.2          last modified: 2019/02/27

if nargin < 5
   Vref.dim(1) = size(img,1);
   Vref.dim(2) = size(img,2);
   Vref.dim(3) = size(img,3);
   Vref.mat    = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1];
end
if nargin < 4
   datatype = 64;
end
if nargin < 3
   description = ' ';
end

% pinfo = Vref.pinfo;

% Vout   = struct('fname', outputfilename,...
%    		      'dim',     Vref.dim(1:3),...
%               'dt',      [datatype 0],...
% 		      'mat',	 Vref.mat,...
%               'pinfo',   pinfo,...
% 		      'descrip', [description]);
if size(img,4) == 1
   Vout   = struct('fname', outputfilename,...
        		   'dim',     Vref.dim(1:3),...
                   'dt',      [datatype 0],...
		           'mat',	 Vref.mat,...
		           'descrip', [description]);
   Vout = spm_write_vol(Vout,img);
else 
   for i = 1:size(img,4)
       Vout(i) = struct('fname', outputfilename,...
                        'dim',     Vref.dim(1:3),...
                        'dt',      [datatype 0],...
                        'mat',	  Vref.mat,...
                        'descrip', [description], ...
                        'n',       [i 1]);
       spm_write_vol(Vout(i),img(:,:,:,i));
   end
end