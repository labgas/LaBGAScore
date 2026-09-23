function outdir = LaBGAScore_export_scaled_contrasts(resultsdir, conidx, outdir, varargin)
% Write the second-level contrast images to disk, one .nii per subject.
%
% :Usage:
% ::
%    outdir = LaBGAScore_export_scaled_contrasts(resultsdir, conidx, outdir, ...)
%
% WHY THIS EXISTS. The decoding wrapper points cfg.files.name at FIRST-LEVEL
% con images, so an SVM run that way analyses raw contrasts while the GLM and
% the PDM analyse contrasts built from z-scored (and, where applicable,
% ComBat-harmonised) condition images. The three analyses then answer the same
% question on different feature spaces. Writing the second-level object out as
% ordinary .nii files lets TDT consume exactly what the GLM used, with the
% provenance visible in the file list rather than injected into passed_data.
%
% :Inputs:
%   **resultsdir**  model results dir holding contrast_data_objects.mat and
%                   image_names_and_setup.mat
%   **conidx**      contrast index into DATA_OBJ_CON*
%   **outdir**      directory to write into (created if absent)
%
% :Optional:
%   **'object'**    which object to export, default 'DATA_OBJ_CONsc'
%                   ('DATA_OBJ_CON' raw, 'DATA_OBJ_CONsc' z-scored conditions,
%                    'DATA_OBJ_CONscc' l2norm-scaled contrasts)
%   **'overwrite'** default false; refuse to clobber an existing non-empty dir
%   **'tag'**       default ''; suffix on the contrast-objects filename, so a model
%                   holding more than one harmonisation path can be exported from the
%                   right one (e.g. 'labelblind' -> contrast_data_objects_labelblind.mat)
%
% :Output:
%   **outdir**      the directory written, one <subjectID>.nii per subject
%
% Subject identity is taken from DAT.imgs, which is the order the object's
% columns are in, so the mapping cannot silently desynchronise.
% -------------------------------------------------------------------------

objname   = 'DATA_OBJ_CONsc';
overwrite = false;
objtag    = '';
for i = 1:2:numel(varargin)
    switch lower(varargin{i})
        case 'object',    objname   = varargin{i+1};
        case 'overwrite', overwrite = varargin{i+1};
        case 'tag',       objtag    = varargin{i+1};
        otherwise, error('unknown option %s', varargin{i});
    end
end

f_obj = fullfile(resultsdir, ['contrast_data_objects' objtag '.mat']);
f_dat = fullfile(resultsdir, 'image_names_and_setup.mat');
if ~exist(f_obj,'file'), error('not found: %s', f_obj); end
if ~exist(f_dat,'file'), error('not found: %s', f_dat); end

L = load(f_obj, objname);
if ~isfield(L, objname), error('%s not in %s', objname, f_obj); end
OBJ = L.(objname);
if conidx > numel(OBJ), error('conidx %d but only %d contrasts', conidx, numel(OBJ)); end
obj = OBJ{conidx};

% Subject identity: DAT.imgs where prep_2 wrote it, else the firstsubjs cell
% that a_set_up_paths saves alongside DAT. proj_discoverie has the first,
% proj_cfs only the second; both are in the object's column order.
S = load(f_dat);
ids = {};
if isfield(S,'DAT') && isfield(S.DAT,'imgs') && numel(S.DAT.imgs) >= conidx
    imgs = S.DAT.imgs{conidx};
    if ischar(imgs), imgs = cellstr(imgs); end
    ids = cell(numel(imgs),1);
    for i = 1:numel(imgs)
        t = regexp(imgs{i}, '(sub-[A-Za-z0-9]+)', 'tokens', 'once');
        if isempty(t), error('no sub-* id in image name: %s', imgs{i}); end
        ids{i} = t{1};
    end
    src_desc = 'DAT.imgs';
elseif isfield(S,'firstsubjs')
    ids = cellstr(string(S.firstsubjs(:)));
    src_desc = 'firstsubjs';
else
    error('neither DAT.imgs nor firstsubjs is available in %s', f_dat);
end
if numel(ids) ~= size(obj.dat,2)
    error('%s has %d entries but the object has %d images', src_desc, numel(ids), size(obj.dat,2));
end
fprintf('subject ids from %s\n', src_desc);
if numel(unique(ids)) ~= numel(ids)
    error('subject ids are not unique - cannot write one file per subject');
end

if ~exist(outdir,'dir')
    mkdir(outdir);
elseif ~overwrite && ~isempty(dir(fullfile(outdir,'*.nii')))
    error('%s already holds .nii files; pass ''overwrite'', true to replace', outdir);
end

fprintf('\nexporting %s contrast %d (%d subjects, %d voxels) to\n  %s\n', ...
    objname, conidx, numel(ids), size(obj.dat,1), outdir);

for i = 1:numel(ids)
    one = get_wh_image(obj, i);
    fn  = fullfile(outdir, [ids{i} '.nii']);
    one.fullpath = fn;
    write(one, 'overwrite');
end

n = numel(dir(fullfile(outdir,'*.nii')));
fprintf('wrote %d files (expected %d)\n', n, numel(ids));
if n ~= numel(ids)
    error('wrote %d files but expected %d', n, numel(ids));
end

% provenance alongside the images
fid = fopen(fullfile(outdir,'EXPORT_INFO.txt'),'w');
fprintf(fid, 'source      : %s\nobject      : %s\ncontrast idx: %d\nsubjects    : %d\nvoxels      : %d\nwritten     : %s\n', ...
    resultsdir, objname, conidx, numel(ids), size(obj.dat,1), datestr(now));
fprintf(fid, '\nOrder follows DAT.imgs, i.e. the object column order.\n');
fclose(fid);

end
