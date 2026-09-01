function PROV = LaBGAScore_prov_snapshot(varargin)
%
%
% *USAGE*
%
% Records exactly which version of every dependency was in place at the
% moment an analysis script ran, and writes that record next to the
% script's results.
%
% The problem it solves: LaBGAS analysis scripts are copied into a study's
% code subdataset and are therefore frozen per project, but their
% dependencies - CanlabCore, CanlabPrivate, the LaBGAS fork of
% CANlab_help_examples, canlab_single_trials and the rest - are shared,
% mutable clones that keep moving. Nothing else in the workflow records
% which commit of those produced a given result.
%
% Call it directly, or let LaBGAScore_prov_publish call it for you as part
% of publishing a script to html.
%
% Two things make the record trustworthy rather than merely reassuring:
%
% # a commit hash alone does not pin a dependency. Repositories here carry
%   uncommitted local edits (CanlabCore currently has five). Those are
%   detected, listed, and their current content is stored so the exact
%   state can be restored later.
% # given a dependency map, the record distinguishes uncommitted edits that
%   this script's call graph actually reaches from ones that are merely
%   present somewhere in the repository. Most of the time the answer is
%   "none of them", which turns an alarming "dirty" into a precise
%   "clean for this script".
%
% Nothing here shells out, so it cannot hang or fail on a broken system()
% and takes a couple of seconds across ~30 repositories.
%
%
% *OPTIONS*
%
% * 'scriptname'    script this snapshot belongs to. Name or full path;
%                   resolved with which() if only a name is given.
% * 'savedir'       where to write the record. Default: the value of
%                   notesdir in the base workspace, else the current dir.
%                   Pass '' to write nothing and only return PROV.
% * 'githubrootdir' default '/data/master_github_repos'
% * 'extraroots'    non-git dependency roots, default {'/opt/KUL_apps/spm12'}
% * 'depmap'        DEP struct from LaBGAScore_dep_map. If omitted and a
%                   scriptname is given, it is computed (see 'usedepmap').
% * 'usedepmap'     compute a dependency map when none is supplied,
%                   default true. Set false for a faster, coarser record.
% * 'checkdirty'    detect uncommitted modifications, default true
% * 'print'         print the record, default true. The printed form is
%                   what publish() captures into the html report.
%
%
% *OUTPUT*
%
% * PROV            struct with fields
%       .deps       table, one row per repository: Repo, Path, Remote,
%                   Branch, Commit, CommitShort, Dirty, NModified,
%                   ModifiedFiles, DirtyRelevant, RelevantModifiedFiles,
%                   UsedByScript, NFilesUsed, ReflogProtected, NReflog
%       .env        MATLAB and platform details
%       .script     name, path, sha1 of the script itself
%       .nonrepo    dependencies that are not git repositories (SPM12)
%       .modified   table of modified dependency files WITH their current
%                   content, so the state can be reconstructed
%       .files      every dependency file the script reaches, with its
%                   repository and current blob sha1
%
%
% *NOTES*
%
% UsedByScript and DirtyRelevant are only meaningful when a dependency map
% is available. Without one they are false throughout and the record falls
% back to repository-level information, which is still correct, just less
% specific.
%
%
% *SEE ALSO*
%
% LaBGAScore_prov_publish, LaBGAScore_prov_gitinfo, LaBGAScore_prov_gitstatus,
% LaBGAScore_dep_map, LaBGAScore_prov_resolve_retrospective
%
% -------------------------------------------------------------------------
%
% modified by: Lukas Van Oudenhove & Claude Opus 5
%
% date:   KU Leuven, August, 2026
%
% -------------------------------------------------------------------------
%
% LaBGAScore_prov_snapshot.m         v1.0
%
% last modified: 2026/08/31
%
%
%% PARSE OPTIONS
% -------------------------------------------------------------------------

p = inputParser;
p.addParameter('scriptname', '', @(x) ischar(x) || isstring(x));
p.addParameter('savedir', '__default__', @(x) ischar(x) || isstring(x));
p.addParameter('githubrootdir', '/data/master_github_repos', @(x) ischar(x) || isstring(x));
p.addParameter('extraroots', {'/opt/KUL_apps/spm12'}, @(x) iscell(x) || isstring(x));
p.addParameter('depmap', [], @(x) isempty(x) || isstruct(x));
p.addParameter('usedepmap', true, @islogical);
p.addParameter('checkdirty', true, @islogical);
p.addParameter('print', true, @islogical);
p.parse(varargin{:});
opt = p.Results;

githubrootdir = char(opt.githubrootdir);


%% RESOLVE THE SCRIPT AND ITS DEPENDENCY MAP
% -------------------------------------------------------------------------

scriptpath = '';
scriptname = char(opt.scriptname);

if ~isempty(scriptname)
    if isfile(scriptname)
        scriptpath = scriptname;
    else
        w = which(scriptname);
        if ~isempty(w) && isfile(w), scriptpath = w; end
    end
    [~, scriptname] = fileparts(scriptname);
end

DEP = opt.depmap;

if isempty(DEP) && opt.usedepmap && ~isempty(scriptpath)
    try
        DEP = LaBGAScore_dep_map(scriptpath, 'print', false);
    catch ME
        warning('LaBGAScore_prov_snapshot:depmap', ...
            'could not build dependency map (%s); recording repository-level provenance only', ...
            ME.message);
        DEP = [];
    end
end

usedfiles = strings(0,1);
mindepthrepo = strings(0,1); mindepthval = zeros(0,1);

if ~isempty(DEP) && isfield(DEP,'files') && ~isempty(DEP.files)
    usedfiles = DEP.files.File;
    % How directly is each repository reached? One called at depth 1 is a real
    % dependency of this script; one appearing only at depth 6 through a chain
    % of ambiguous names is incidental. The dirty check below still uses the
    % full closure - over-approximating is the safe direction there - but the
    % report needs to tell the two apart.
    e = DEP.edges(DEP.edges.Repo ~= "",:);
    if ~isempty(e)
        [mindepthrepo, ~, ix] = unique(e.Repo);
        mindepthval = accumarray(ix, e.Depth, [], @min);
    end
end


%% SNAPSHOT EACH REPOSITORY
% -------------------------------------------------------------------------

d = dir(githubrootdir);
d = d([d.isdir]);
d = d(~startsWith({d.name},'.'));

repo = strings(0,1); repopath = strings(0,1); remote = strings(0,1);
branch = strings(0,1); commit = strings(0,1); commitshort = strings(0,1);
dirty = false(0,1); nmodified = zeros(0,1); modifiedfiles = strings(0,1);
dirtyrelevant = false(0,1); relevantmodified = strings(0,1);
usedbyscript = false(0,1); nfilesused = zeros(0,1); mindepth = zeros(0,1);
reflogprotected = false(0,1); nreflog = zeros(0,1);

modrepo = strings(0,1); modpath = strings(0,1); modindexsha = strings(0,1);
modworksha = strings(0,1); modcontent = strings(0,1); modrelevant = false(0,1);

for r = 1:numel(d)

    thisrepo = d(r).name;
    thispath = fullfile(githubrootdir, thisrepo);

    G = LaBGAScore_prov_gitinfo(thispath);
    if ~G.isrepo, continue, end

    % which of this script's dependency files live in this repo
    inrepo = usedfiles(startsWith(usedfiles, string(thispath) + filesep));
    relused = erase(inrepo, string(thispath) + filesep);

    % dirty state
    thismod = strings(0,1);
    if opt.checkdirty
        S = LaBGAScore_prov_gitstatus(thispath);
        if S.ok
            thismod = [S.modified; S.deleted];
        end
    end

    thisrelevant = intersect(thismod, relused);

    repo(end+1,1) = string(thisrepo); %#ok<AGROW>
    repopath(end+1,1) = string(thispath); %#ok<AGROW>
    remote(end+1,1) = string(G.remote); %#ok<AGROW>
    branch(end+1,1) = string(G.branch); %#ok<AGROW>
    commit(end+1,1) = string(G.commit); %#ok<AGROW>
    commitshort(end+1,1) = string(G.commit_short); %#ok<AGROW>
    dirty(end+1,1) = ~isempty(thismod); %#ok<AGROW>
    nmodified(end+1,1) = numel(thismod); %#ok<AGROW>
    modifiedfiles(end+1,1) = strjoin(thismod, '; '); %#ok<AGROW>
    dirtyrelevant(end+1,1) = ~isempty(thisrelevant); %#ok<AGROW>
    relevantmodified(end+1,1) = strjoin(thisrelevant, '; '); %#ok<AGROW>
    usedbyscript(end+1,1) = ~isempty(inrepo); %#ok<AGROW>
    nfilesused(end+1,1) = numel(inrepo); %#ok<AGROW>
    md = mindepthval(mindepthrepo == string(thisrepo));
    if isempty(md), mindepth(end+1,1) = NaN; else, mindepth(end+1,1) = md(1); end %#ok<AGROW>
    reflogprotected(end+1,1) = G.reflog_protected; %#ok<AGROW>
    nreflog(end+1,1) = height(G.reflog); %#ok<AGROW>

    % keep the content of modified files that the script actually reaches,
    % so the exact state can be restored even though it was never committed
    for k = 1:numel(thisrelevant)
        fp = fullfile(thispath, char(thisrelevant(k)));
        modrepo(end+1,1) = string(thisrepo); %#ok<AGROW>
        modpath(end+1,1) = thisrelevant(k); %#ok<AGROW>
        modindexsha(end+1,1) = ""; %#ok<AGROW>
        modrelevant(end+1,1) = true; %#ok<AGROW>
        if isfile(fp)
            modworksha(end+1,1) = string(local_filesha(fp)); %#ok<AGROW>
            modcontent(end+1,1) = string(fileread(fp)); %#ok<AGROW>
        else
            modworksha(end+1,1) = ""; %#ok<AGROW>
            modcontent(end+1,1) = "<deleted>"; %#ok<AGROW>
        end
    end

end

PROV.deps = table(repo, repopath, remote, branch, commit, commitshort, ...
    dirty, nmodified, modifiedfiles, dirtyrelevant, relevantmodified, ...
    usedbyscript, nfilesused, mindepth, reflogprotected, nreflog, ...
    'VariableNames', {'Repo','Path','Remote','Branch','Commit','CommitShort', ...
        'Dirty','NModified','ModifiedFiles','DirtyRelevant','RelevantModifiedFiles', ...
        'UsedByScript','NFilesUsed','MinDepth','ReflogProtected','NReflog'});

PROV.deps = sortrows(PROV.deps, {'MinDepth','NFilesUsed'}, {'ascend','descend'});

PROV.modified = table(modrepo, modpath, modindexsha, modworksha, modrelevant, modcontent, ...
    'VariableNames', {'Repo','Path','IndexSHA','WorkingSHA','Relevant','Content'});


%% RECORD THE DEPENDENCY FILES THEMSELVES
% -------------------------------------------------------------------------

if isempty(usedfiles)
    PROV.files = table(strings(0,1), strings(0,1), strings(0,1), ...
        'VariableNames', {'File','Repo','SHA'});
else
    sha = strings(numel(usedfiles),1);
    for k = 1:numel(usedfiles)
        if isfile(usedfiles(k))
            sha(k) = string(local_filesha(char(usedfiles(k))));
        end
    end
    PROV.files = table(usedfiles, DEP.files.Repo, sha, ...
        'VariableNames', {'File','Repo','SHA'});
end


%% ENVIRONMENT AND SCRIPT
% -------------------------------------------------------------------------

PROV.env = struct();
PROV.env.timestamp = string(datetime('now','TimeZone','local', ...
    'Format','yyyy-MM-dd''T''HH:mm:ssXXX'));
PROV.env.matlab_version = string(version);
PROV.env.matlab_release = string(version('-release'));
PROV.env.computer = string(computer);
PROV.env.hostname = string(local_hostname);
PROV.env.user = string(getenv('USER'));
PROV.env.pwd = string(pwd);

% Screen geometry, because figure captures depend on it. These scripts size
% figures for publish() either by maximising the window or via
% plugin_set_figure_size, and BOTH are bounded by the display of whichever
% client was connected when the script ran: publish captures what is on
% screen, so a figure larger than the display comes back at display size,
% with a different aspect ratio than intended. Recording the screen makes
% that difference visible across machines instead of merely puzzling.
try
    ss = get(0,'ScreenSize');
    PROV.env.screen_size = string(sprintf('%dx%d', ss(3), ss(4)));
    PROV.env.screen_dpi = get(0,'ScreenPixelsPerInch');
    PROV.env.max_figure_inches = string(sprintf('%.1fx%.1f', ...
        ss(3)/PROV.env.screen_dpi, ss(4)/PROV.env.screen_dpi));
catch
    PROV.env.screen_size = "unknown";
    PROV.env.screen_dpi = NaN;
    PROV.env.max_figure_inches = "unknown";
end

PROV.script = struct();
PROV.script.name = string(scriptname);
PROV.script.path = string(scriptpath);
if ~isempty(scriptpath) && isfile(scriptpath)
    PROV.script.sha = string(local_filesha(scriptpath));
else
    PROV.script.sha = "";
end


%% NON-GIT DEPENDENCIES
% -------------------------------------------------------------------------

extraroots = cellstr(opt.extraroots);
nrname = strings(0,1); nrpath = strings(0,1); nrversion = strings(0,1);

for k = 1:numel(extraroots)
    thispath = extraroots{k};
    if ~isfolder(thispath), continue, end
    [~, thisname] = fileparts(thispath);
    nrname(end+1,1) = string(thisname); %#ok<AGROW>
    nrpath(end+1,1) = string(thispath); %#ok<AGROW>
    nrversion(end+1,1) = string(local_contents_version(thispath)); %#ok<AGROW>
end

PROV.nonrepo = table(nrname, nrpath, nrversion, ...
    'VariableNames', {'Name','Path','Version'});


%% WRITE THE RECORD
% -------------------------------------------------------------------------

savedir = char(opt.savedir);

if strcmp(savedir,'__default__')
    savedir = '';
    try
        if evalin('base','exist(''notesdir'',''var'')')
            savedir = evalin('base','notesdir');
        end
    catch
    end
    if isempty(savedir), savedir = pwd; end
end

PROV.savedir = string(savedir);
PROV.files_written = strings(0,1);

if ~isempty(savedir)

    if ~isfolder(savedir), mkdir(savedir); end

    stamp = char(datetime('now','Format','yyyyMMdd''T''HHmmss'));
    if isempty(scriptname), base = 'provenance'; else, base = scriptname; end
    stem = fullfile(savedir, sprintf('%s_provenance_%s', base, stamp));

    save([stem '.mat'], 'PROV', '-v7.3');
    PROV.files_written(end+1,1) = string([stem '.mat']);

    % small, text, and well under the 100kb git-annex threshold used by the
    % project .gitattributes, so it stays in git proper and stays readable
    % without a datalad get
    T = PROV.deps;
    writetable(T, [stem '.tsv'], 'FileType','text', 'Delimiter','\t');
    PROV.files_written(end+1,1) = string([stem '.tsv']);

    if ~isempty(PROV.modified)
        writetable(removevars(PROV.modified,'Content'), [stem '_modified.tsv'], ...
            'FileType','text', 'Delimiter','\t');
        PROV.files_written(end+1,1) = string([stem '_modified.tsv']);
    end

end


%% PRINT THE RECORD
% -------------------------------------------------------------------------

if opt.print
    local_print(PROV, DEP);
end

end


%% LOCAL FUNCTIONS
% -------------------------------------------------------------------------

function local_print(PROV, DEP)

dashes = repmat('-',1,75);

fprintf('\n%s\n', dashes);
fprintf('PROVENANCE\n');
fprintf('%s\n', dashes);

fprintf('script    %s\n', PROV.script.name);
fprintf('run       %s on %s by %s\n', PROV.env.timestamp, PROV.env.hostname, PROV.env.user);
fprintf('MATLAB    %s (R%s), %s\n', PROV.env.matlab_version, PROV.env.matlab_release, PROV.env.computer);
if isfield(PROV.env,'screen_size')
    fprintf('screen    %s @ %g DPI (fits a figure up to %s inches)\n', ...
        PROV.env.screen_size, PROV.env.screen_dpi, PROV.env.max_figure_inches);
end

for k = 1:height(PROV.nonrepo)
    fprintf('%-9s %s\n', PROV.nonrepo.Name(k), PROV.nonrepo.Version(k));
end

used = PROV.deps(PROV.deps.UsedByScript,:);

if isempty(DEP) || isempty(used)
    fprintf('\nNo dependency map available - reporting all repositories.\n');
    used = PROV.deps;
    showall = true;
    nindirect = 0;
else
    showall = false;
    % direct callees, plus anything with a relevant uncommitted change, lead;
    % repos reachable only through long ambiguous chains are counted not listed
    isprimary = used.MinDepth <= 2 | used.DirtyRelevant;
    nindirect = sum(~isprimary);
    used = used(isprimary,:);
    fprintf('\nRepositories this script reaches (%d of %d):\n\n', height(used), height(PROV.deps));
end

if showall
    fprintf('  %-28s %-10s %-9s  %s\n', 'REPOSITORY','COMMIT','BRANCH','STATE');
else
    fprintf('  %-28s %-10s %-9s %5s %6s  %s\n', 'REPOSITORY','COMMIT','BRANCH','DEPTH','FILES','STATE');
end
for k = 1:height(used)
    if used.DirtyRelevant(k)
        state = sprintf('%d MODIFIED file(s) THIS SCRIPT USES', ...
            numel(split(used.RelevantModifiedFiles(k),'; ')));
    elseif used.Dirty(k) && showall
        % without a dependency map there is no basis for saying the changes
        % do not matter here, so do not imply it
        state = sprintf('%d modified file(s); relevance unknown (no dependency map)', ...
            used.NModified(k));
    elseif used.Dirty(k)
        state = sprintf('clean for this script (%d modified elsewhere)', used.NModified(k));
    else
        state = 'clean';
    end
    if showall
        fprintf('  %-28s %-10s %-9s  %s\n', used.Repo(k), used.CommitShort(k), ...
            used.Branch(k), state);
    else
        fprintf('  %-28s %-10s %-9s %5d %6d  %s\n', used.Repo(k), used.CommitShort(k), ...
            used.Branch(k), used.MinDepth(k), used.NFilesUsed(k), state);
    end
end

if nindirect > 0
    fprintf('  (%d further repositories reachable only indirectly, at depth >2;\n', nindirect);
    fprintf('   full list in PROV.deps)\n');
end

flagged = used(used.DirtyRelevant,:);
if ~isempty(flagged)
    fprintf('\n  *** UNCOMMITTED CHANGES IN FILES THIS SCRIPT USES ***\n');
    fprintf('  The commit hashes above do NOT fully describe what ran.\n');
    for k = 1:height(flagged)
        files = split(flagged.RelevantModifiedFiles(k),'; ');
        for i = 1:numel(files)
            fprintf('    %s: %s\n', flagged.Repo(k), files(i));
        end
    end
    fprintf('  Current content of these files is stored in the .mat record.\n');
end

if ~showall
    unprot = PROV.deps(PROV.deps.UsedByScript & ~PROV.deps.ReflogProtected,:);
    if ~isempty(unprot)
        fprintf('\n  Note: %d of these repositories can still lose their reflog to git gc,\n', height(unprot));
        fprintf('  which would make this run unreconstructable later. Run\n');
        fprintf('  labgascore_prov_protect_reflogs.sh to fix.\n');
    end
end

if ~isempty(PROV.files_written)
    fprintf('\nRecord written to:\n');
    for k = 1:numel(PROV.files_written)
        fprintf('  %s\n', PROV.files_written(k));
    end
end

fprintf('%s\n\n', dashes);

end


function h = local_filesha(fp)
% git-compatible blob sha1, so it can be compared against git object ids

fid = fopen(fp,'r');
if fid < 0, h = ''; return, end
data = fread(fid, Inf, '*uint8');
fclose(fid);

md = java.security.MessageDigest.getInstance('SHA-1');
md.update([uint8(sprintf('blob %d', numel(data))), uint8(0)]);
md.update(data);

h = char(lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2)', 1, [])));

end


function h = local_hostname
try
    h = char(java.net.InetAddress.getLocalHost().getHostName());
catch
    h = getenv('HOSTNAME');
end
if isempty(h), h = 'unknown'; end
end


function v = local_contents_version(rootpath)
% version string from a toolbox Contents.m, e.g. SPM12's "Version 7771"

v = 'unknown';
cf = fullfile(rootpath,'Contents.m');
if ~isfile(cf), return, end

lines = splitlines(string(fileread(cf)));
for k = 1:min(20,numel(lines))
    t = regexp(lines(k), '^%\s*(Version\s+.+)$', 'tokens', 'once');
    if ~isempty(t), v = strtrim(char(t)); return, end
end

end
