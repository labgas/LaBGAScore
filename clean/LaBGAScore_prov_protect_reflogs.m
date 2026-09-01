function report = LaBGAScore_prov_protect_reflogs(githubrootdir, archivedir)
%
%
% *USAGE*
%
% Archives the git reflog of every repository cloned under githubrootdir,
% and reports whether each one is protected against reflog expiry.
%
% The reflog is what lets LaBGAScore_prov_resolve_retrospective work out
% which commit of a dependency was checked out when a past analysis was
% run. It is the ONLY record of that: a commit's author/committer date says
% nothing about when your clone actually moved to it. Git's default
% gc.reflogExpire is 90 days and expiry happens silently during garbage
% collection, so once pruned the provenance of past analyses is gone.
%
% RUN THIS PERIODICALLY. The archived .tsv files, not the live reflogs, are
% what the retrospective resolver reads, so the record survives a later
% re-clone or a gc.
%
% This function does NOT shell out - it reads .git/logs/HEAD and .git/config
% directly via LaBGAScore_prov_gitinfo. It therefore also cannot SET the
% expiry config. Any repository found unprotected is listed at the end
% along with the exact shell command to fix it; there is also a ready-made
% script next to this file:
%
%   bash labgascore_prov_protect_reflogs.sh [githubrootdir]
%
% Nothing is ever deleted or overwritten: each run writes a new timestamped
% archive alongside the existing ones.
%
%
% *OPTIONS*
%
% * githubrootdir   directory holding the local clones, default
%                   '/data/master_github_repos'
% * archivedir      where to write the archives, default
%                   <githubrootdir>/.labgas_provenance
%
%
% *OUTPUT*
%
% * report          table, one row per directory found, with columns
%                   Repo, Path, IsGitRepo, Protected, Branch, Commit,
%                   Remote, NEntries, FirstEntry, LastEntry, Archive, Status
%
%
% *SEE ALSO*
%
% LaBGAScore_prov_gitinfo, LaBGAScore_prov_resolve_retrospective
%
% -------------------------------------------------------------------------
%
% modified by: Lukas Van Oudenhove & Claude Opus 5
%
% date:   KU Leuven, August, 2026
%
% -------------------------------------------------------------------------
%
% LaBGAScore_prov_protect_reflogs.m         v1.1
%
% last modified: 2026/08/31
%
%
%% SET DEFAULTS
% -------------------------------------------------------------------------

if nargin < 1 || isempty(githubrootdir)
    githubrootdir = '/data/master_github_repos';
end

if nargin < 2 || isempty(archivedir)
    archivedir = fullfile(githubrootdir,'.labgas_provenance');
end

if ~isfolder(githubrootdir)
    error('githubrootdir %s does not exist', githubrootdir);
end

if ~isfolder(archivedir)
    mkdir(archivedir);
    fprintf('\nCreated reflog archive directory %s\n', archivedir);
end

stamp = char(datetime('now','Format','yyyyMMdd''T''HHmmss'));


%% SCAN AND ARCHIVE
% -------------------------------------------------------------------------

d = dir(githubrootdir);
d = d([d.isdir]);
d = d(~startsWith({d.name},'.'));           % skips . .. and .labgas_provenance

repo = strings(0,1); repopath = strings(0,1); isgitrepo = false(0,1);
protected = false(0,1); branch = strings(0,1); commit = strings(0,1);
remote = strings(0,1); nentries = zeros(0,1); firstentry = strings(0,1);
lastentry = strings(0,1); archive = strings(0,1); status = strings(0,1);

fprintf('\nArchiving reflogs for repositories under %s\n\n', githubrootdir);

for r = 1:numel(d)

    thisrepo = d(r).name;
    thispath = fullfile(githubrootdir, thisrepo);

    G = LaBGAScore_prov_gitinfo(thispath);

    repo(end+1,1) = string(thisrepo); %#ok<AGROW>
    repopath(end+1,1) = string(thispath); %#ok<AGROW>
    isgitrepo(end+1,1) = G.isrepo; %#ok<AGROW>
    protected(end+1,1) = G.reflog_protected; %#ok<AGROW>
    branch(end+1,1) = string(G.branch); %#ok<AGROW>
    commit(end+1,1) = string(G.commit_short); %#ok<AGROW>
    remote(end+1,1) = string(G.remote); %#ok<AGROW>

    if ~G.isrepo
        nentries(end+1,1) = NaN; firstentry(end+1,1) = ""; %#ok<AGROW>
        lastentry(end+1,1) = ""; archive(end+1,1) = ""; %#ok<AGROW>
        status(end+1,1) = string(G.error); %#ok<AGROW>
        fprintf('  %-32s skipped (%s)\n', thisrepo, G.error);
        continue
    end

    n = height(G.reflog);
    nentries(end+1,1) = n; %#ok<AGROW>

    if n == 0
        firstentry(end+1,1) = ""; lastentry(end+1,1) = ""; %#ok<AGROW>
        archive(end+1,1) = ""; %#ok<AGROW>
        status(end+1,1) = "empty reflog"; %#ok<AGROW>
        fprintf('  %-32s empty reflog - past checkouts NOT reconstructable\n', thisrepo);
        continue
    end

    archivefile = fullfile(archivedir, sprintf('reflog_%s_%s.tsv', thisrepo, stamp));

    T = G.reflog;
    T.Time = string(T.Time, 'yyyy-MM-dd''T''HH:mm:ss');      % local, stable ISO 8601
    writetable(T, archivefile, 'FileType','text', 'Delimiter','\t');

    firstentry(end+1,1) = T.Time(1); %#ok<AGROW>
    lastentry(end+1,1) = T.Time(end); %#ok<AGROW>
    archive(end+1,1) = string(archivefile); %#ok<AGROW>

    if G.reflog_protected
        status(end+1,1) = "archived, expiry disabled"; %#ok<AGROW>
    else
        status(end+1,1) = "archived, EXPIRY STILL ACTIVE"; %#ok<AGROW>
    end

    fprintf('  %-32s %4d entries, %s .. %s%s\n', thisrepo, n, ...
        extractBefore(T.Time(1),11), extractBefore(T.Time(end),11), ...
        local_flag(G.reflog_protected));

end

report = table(repo, repopath, isgitrepo, protected, branch, commit, remote, ...
    nentries, firstentry, lastentry, archive, status, ...
    'VariableNames', {'Repo','Path','IsGitRepo','Protected','Branch','Commit', ...
                      'Remote','NEntries','FirstEntry','LastEntry','Archive','Status'});


%% PRINT REPORT
% -------------------------------------------------------------------------

ngit = sum(report.IsGitRepo);
narch = sum(report.NEntries > 0);

fprintf('\n%d directories scanned, %d git repositories, %d reflogs archived\n', ...
    height(report), ngit, narch);
fprintf('Archives in %s\n', archivedir);

unprot = report(report.IsGitRepo & ~report.Protected, :);

if isempty(unprot)
    fprintf('\nAll %d repositories have reflog expiry disabled.\n\n', ngit);
else
    fprintf('\n*** %d repository(ies) can still LOSE their reflog to git gc ***\n', height(unprot));
    for k = 1:height(unprot)
        fprintf('  %s\n', unprot.Repo(k));
    end
    fprintf('\nFix from a shell (this function deliberately does not shell out):\n');
    fprintf('  bash %s %s\n\n', ...
        fullfile(fileparts(mfilename('fullpath')),'labgascore_prov_protect_reflogs.sh'), ...
        githubrootdir);
end

empty = report(report.IsGitRepo & report.NEntries == 0, :);
if ~isempty(empty)
    fprintf('Note: %d repository(ies) have an empty reflog. Runs against those\n', height(empty));
    fprintf('resolve as committer-date_APPROX rather than exactly:\n');
    for k = 1:height(empty)
        fprintf('  %s\n', empty.Repo(k));
    end
    fprintf('\n');
end

end


%% LOCAL FUNCTIONS
% -------------------------------------------------------------------------

function s = local_flag(isprotected)
if isprotected
    s = '';
else
    s = '   <-- EXPIRY STILL ACTIVE';
end
end
