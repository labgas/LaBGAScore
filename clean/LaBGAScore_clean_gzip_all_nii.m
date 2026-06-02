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
%   if not specified or left empty, all .nii files will be gzipped
% 3. delete unzipped (true/false)
%   specify whether or not to delete the unzipped files
%   if not specified, default = false
%   false is useful if you want to use datalad remove to delete all
%   unzipped files in one go after running this function - recommended
%
% OUTPUT
% cell array with paths of zipped .nii files
%
% AUTHOR
% Lukas Van Oudenhove, KU Leuven, May 2026
%
% ..

%% PARSE VARARGIN

if isempty(varargin)
    input_folder = pwd;
else
    input_folder = varargin{1};
end

if length(varargin) > 1 && ~isempty(varargin{2})
    filter_label = varargin{2};
end

do_delete = false;

if length(varargin) > 2
    do_delete = varargin{3};
end

%% CORE FUNCTION

cd(input_folder);

if ~exist("filter_label","var")
    niidir = dir('**/*.nii');
else
    niidir = dir(['**/*' filter_label '*.nii']);
end

niis2gzip = cell(size(niidir,1),1);

LaBGAScore_smart_parallel_pool_setup

if do_delete

    parfor i = 1:size(niidir,1)
        niis2gzip{i} = fullfile(niidir(i).folder,niidir(i).name);
        gzip(niis2gzip{i});
        delete(niis2gzip{i});
    end

else 
    
    parfor i = 1:size(niidir,1)
        niis2gzip{i} = fullfile(niidir(i).folder,niidir(i).name);
        gzip(niis2gzip{i});
    end
    
end

%% OUTPUT

gzippedniis = niis2gzip;

fprintf('\ngzipped and deleted %s .nii files in %s\n\n', num2str(size(gzippedniis,1)), input_folder);

end