function LaBGAScore_clean_sourcedata(path_to_superdataset)
%
% cleans sourcedata subdataset depending on whether DICOMs are gitignored
% or not
%
% see 
%
% OPTIONAL INPUT
% path_to_superdataset : absolute or relative path to superdataset
% if not specified, present working directory is used
%
% AUTHOR
% Lukas Van Oudenhove, KU Leuven, Jan 2026
%
% ..

if exist('path_to_superdataset','var')
    sourcedir = fullfile(path_to_superdataset,'sourcedata');
else
    sourcedir = fullfile(pwd,'sourcedata');
end

cd(sourcedir);

[~,cmdout] = system('cat .gitignore');

    if contains(cmdout,'**/DICOM') % DICOMS are gitignored -> simply delete, but keep folder structure
       !find -type f -delete
    else                           % DICOMS are not gitignored -> drop annexed content (but classic DICOMS stored in git rather than git annex will remain present I guess)
       !datalad drop --reckless availability 
    end

end