function gzippedniis = LaBGAScore_clean_gzip_all_nii(varargin)
%
% gzips .nii files in input folder recursively
% optional specification of label to limit zipping to for example beta.nii
% files
%
% INPUTS
% 1. input folder
%   specify absolute path to folder in which you want to gzip recursively
%   if not specified, pwd will be used
% 2. filter
%   specify label to filter, e.g. 'beta', to only gzip beta images
%   if not specified, all .nii files will be gzipped
%
% OUTPUT
% cell array with paths of zipped .nii files
%
% AUTHOR
% Lukas Van Oudenhove, KU Leuven, Feb 2023
%
% ..

%% PARSE VARARGIN

if isempty(varargin)
    input_folder = pwd;
else
    input_folder = varargin{1};
end

if length(varargin) > 1
    filter_label = varargin{2};
end

%% CORE FUNCTION
cd(input_folder);
if ~exist("filter_label","var")
    niidir = dir('**/*.nii');
else
    niidir = dir(['**/*' filter_label '*.nii']);
end

niis2gzip = cell(size(niidir,1),1);

parfor i = 1:size(niidir,1)
    niis2gzip{i} = fullfile(niidir(i).folder,niidir(i).name);
    gzip(niis2gzip{i});
    delete(niis2gzip{i});
end

%% OUTPUT
gzippedniis = niis2gzip;

fprintf('\ngzipped and deleted %s .nii files in %s\n\n', num2str(size(gzippedniis,1)), input_folder);

end