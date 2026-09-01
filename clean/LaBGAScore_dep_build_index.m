function IDX = LaBGAScore_dep_build_index(varargin)
%
%
% *USAGE*
%
% Builds an index of every callable file (.m, .mex*, .p) under the given
% dependency roots, mapping each callable NAME to ALL of its candidate
% definitions. LaBGAScore_dep_map resolves identifiers against this index.
%
% The index exists because which() cannot be trusted for this job. which()
% returns exactly one answer, determined by the current MATLAB path order,
% and on this setup it is demonstrably the wrong one for the calls that
% matter most:
%
%   which('predict')   -> a liblinear .mexa64 in SPM's decoding toolbox
%   which('ttest')     -> MATLAB's own Statistics Toolbox
%   which('threshold') -> CanlabCore/@image_vector/threshold.m
%
% whereas the second-level scripts actually call @fmri_data/predict,
% @fmri_data/ttest and @statistic_image/threshold. CanlabCore alone has 15
% @class directories. Enumerating every candidate and reporting ambiguity
% is therefore the only honest resolution strategy.
%
% The index is cached, keyed on the roots, and only rebuilt when asked.
%
%
% *OPTIONS*
%
% * 'roots'         cellstr of directories to index. Default: every git
%                   repository directly under githubrootdir, plus SPM12,
%                   plus the MATLAB toolbox root.
% * 'githubrootdir' default '/data/master_github_repos'
% * 'extraroots'    cellstr appended to the default roots, default
%                   {'/opt/KUL_apps/spm12'}
% * 'cachefile'     default <clean dir>/.dep_index_<hostname>.mat
% * 'rebuild'       logical, force a rebuild, default false
% * 'print'         logical, default true
%
%
% *OUTPUT*
%
% * IDX             struct with fields
%                   .files  table, one row per callable file, columns
%                           Name, File, Root, Repo, Class, Kind, IsPrivate
%                   .map    containers.Map from lowercase name to the row
%                           indices of .files that define it
%                   .roots  the roots indexed
%                   .built  datetime the index was built
%
%
% *NOTES*
%
% Kind is determined from the first non-comment line: 'classdef',
% 'function' or 'script'. A file inside an @class directory is 'method'
% regardless. Files inside a private/ directory are flagged IsPrivate,
% because they are only callable from their sibling directory.
%
%
% *SEE ALSO*
%
% LaBGAScore_dep_map, LaBGAScore_prov_snapshot
%
% -------------------------------------------------------------------------
%
% modified by: Lukas Van Oudenhove & Claude Opus 5
%
% date:   KU Leuven, August, 2026
%
% -------------------------------------------------------------------------
%
% LaBGAScore_dep_build_index.m         v1.0
%
% last modified: 2026/08/31
%
%
%% PARSE OPTIONS
% -------------------------------------------------------------------------

p = inputParser;
p.addParameter('roots', {}, @(x) iscell(x) || ischar(x) || isstring(x));
p.addParameter('githubrootdir', '/data/master_github_repos', @(x) ischar(x) || isstring(x));
p.addParameter('extraroots', {'/opt/KUL_apps/spm12'}, @(x) iscell(x) || ischar(x) || isstring(x));
p.addParameter('cachefile', '', @(x) ischar(x) || isstring(x));
p.addParameter('rebuild', false, @islogical);
p.addParameter('print', true, @islogical);
p.parse(varargin{:});
opt = p.Results;

githubrootdir = char(opt.githubrootdir);

if isempty(opt.roots)
    roots = local_default_roots(githubrootdir, cellstr(opt.extraroots));
else
    roots = cellstr(opt.roots);
end

roots = roots(cellfun(@isfolder, roots));

if isempty(opt.cachefile)
    try
        host = char(java.net.InetAddress.getLocalHost().getHostName());
    catch
        host = getenv('HOSTNAME');
    end
    if isempty(host), host = 'unknown'; end
    host = matlab.lang.makeValidName(strtok(host,'.'));
    cachefile = fullfile(fileparts(mfilename('fullpath')), ['.dep_index_' host '.mat']);
else
    cachefile = char(opt.cachefile);
end


%% RETURN CACHED INDEX IF IT MATCHES
% -------------------------------------------------------------------------

if ~opt.rebuild && isfile(cachefile)
    S = load(cachefile,'IDX');
    if isequal(sort(S.IDX.roots(:)), sort(roots(:)))
        IDX = S.IDX;
        if opt.print
            fprintf('\nLoaded cached dependency index (%d files, %d names, built %s)\n', ...
                height(IDX.files), IDX.map.Count, char(IDX.built));
            fprintf('Pass ''rebuild'',true to rebuild it.\n\n');
        end
        return
    end
end


%% SCAN THE ROOTS
% -------------------------------------------------------------------------

if opt.print
    fprintf('\nBuilding dependency index over %d roots (this takes a minute)\n', numel(roots));
end

name = strings(0,1); file = strings(0,1); rootcol = strings(0,1);
repo = strings(0,1); class = strings(0,1); kind = strings(0,1);
isprivate = false(0,1);

for r = 1:numel(roots)

    thisroot = roots{r};
    ismatlabroot = startsWith(thisroot, matlabroot);

    d = [dir(fullfile(thisroot,'**','*.m')); ...
         dir(fullfile(thisroot,'**','*.p')); ...
         dir(fullfile(thisroot,'**','*.mex*'))];
    d = d(~[d.isdir]);

    if opt.print
        fprintf('  %-46s %6d files\n', local_shorten(thisroot), numel(d));
    end

    for k = 1:numel(d)

        fullpath = fullfile(d(k).folder, d(k).name);
        [~, basename, ext] = fileparts(d(k).name);

        % Contents.m is documentation, never a callable
        if strcmpi(basename,'Contents'), continue, end

        folderparts = strsplit(d(k).folder, filesep);

        % @class directory -> a method of that class
        atdirs = folderparts(startsWith(folderparts,'@'));
        if ~isempty(atdirs)
            thisclass = extractAfter(atdirs{end},1);
        else
            thisclass = '';
        end

        thisprivate = any(strcmp(folderparts,'private'));

        if strcmpi(ext,'.p')
            thiskind = 'pcode';
        elseif startsWith(ext,'.mex')
            thiskind = 'mex';
        elseif ~isempty(thisclass)
            thiskind = 'method';
        elseif ismatlabroot
            % MathWorks' own code is only ever needed to recognise a name as
            % a builtin so it can be filtered out; reading ~40k files just to
            % label them function vs script would dominate the build
            thiskind = 'toolbox';
        else
            thiskind = local_kind_from_source(fullpath);
        end

        % the repo/toolbox this file belongs to: first component below the root
        rel = erase(fullpath,[thisroot filesep]);
        relparts = strsplit(rel, filesep);
        if strcmp(thisroot, githubrootdir) && numel(relparts) > 1
            thisrepo = relparts{1};
        else
            [~, thisrepo] = fileparts(thisroot);
        end

        name(end+1,1) = string(basename); %#ok<AGROW>
        file(end+1,1) = string(fullpath); %#ok<AGROW>
        rootcol(end+1,1) = string(thisroot); %#ok<AGROW>
        repo(end+1,1) = string(thisrepo); %#ok<AGROW>
        class(end+1,1) = string(thisclass); %#ok<AGROW>
        kind(end+1,1) = string(thiskind); %#ok<AGROW>
        isprivate(end+1,1) = thisprivate; %#ok<AGROW>

    end

end

IDX.files = table(name, file, rootcol, repo, class, kind, isprivate, ...
    'VariableNames', {'Name','File','Root','Repo','Class','Kind','IsPrivate'});
IDX.roots = roots;
IDX.built = datetime('now');


%% BUILD THE NAME -> ROWS MAP
% -------------------------------------------------------------------------

lowernames = lower(IDX.files.Name);
[uniquenames, ~, ix] = unique(lowernames);

IDX.map = containers.Map('KeyType','char','ValueType','any');
for k = 1:numel(uniquenames)
    IDX.map(char(uniquenames(k))) = find(ix == k);
end

save(cachefile,'IDX','-v7.3');


%% PRINT REPORT
% -------------------------------------------------------------------------

if opt.print
    fprintf('\n%d callable files indexed, %d distinct names\n', height(IDX.files), IDX.map.Count);
    nmethod = sum(IDX.files.Kind == "method");
    nambig = 0;
    k = keys(IDX.map);
    for i = 1:numel(k)
        if numel(IDX.map(k{i})) > 1, nambig = nambig + 1; end
    end
    fprintf('%d class methods, %d names with more than one definition\n', nmethod, nambig);
    fprintf('Cached to %s\n\n', cachefile);
end

end


%% LOCAL FUNCTIONS
% -------------------------------------------------------------------------

function roots = local_default_roots(githubrootdir, extraroots)
% every git repo directly under githubrootdir, plus the extras, plus MATLAB

roots = {};

d = dir(githubrootdir);
d = d([d.isdir]);
d = d(~startsWith({d.name},'.'));

for k = 1:numel(d)
    thispath = fullfile(githubrootdir, d(k).name);
    if isfolder(fullfile(thispath,'.git')) || isfile(fullfile(thispath,'.git'))
        roots{end+1} = thispath; %#ok<AGROW>
    end
end

roots = [roots, extraroots(:)', {fullfile(matlabroot,'toolbox')}];

end


function kind = local_kind_from_source(fullpath)
% first non-comment, non-blank line decides function vs classdef vs script

kind = 'script';

fid = fopen(fullpath,'r');
if fid < 0, return, end
c = onCleanup(@() fclose(fid));

for k = 1:200
    tline = fgetl(fid);
    if ~ischar(tline), return, end
    tline = strtrim(tline);
    if isempty(tline) || startsWith(tline,'%'), continue, end
    if ~isempty(regexp(tline,'^classdef\>','once'))
        kind = 'classdef';
    elseif ~isempty(regexp(tline,'^function\>','once'))
        kind = 'function';
    end
    return
end

end


function s = local_shorten(p)
if strlength(p) > 46
    s = ['...' char(extractAfter(string(p), strlength(p)-43))];
else
    s = p;
end
end
