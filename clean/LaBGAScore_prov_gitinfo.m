function G = LaBGAScore_prov_gitinfo(repopath)
%
%
% *USAGE*
%
% Reads the state of a git repository - current commit, branch, remote,
% full reflog and reflog-expiry settings - by parsing the plain-text files
% under .git directly, WITHOUT shelling out to git.
%
% Everything this function returns lives in documented plain-text files:
%
%   .git/HEAD              current branch ref, or a bare sha if detached
%   .git/refs/heads/<b>    the commit that branch points at
%   .git/packed-refs       same, for refs that have been packed away
%   .git/logs/HEAD         the reflog, one line per HEAD movement
%   .git/config            remotes, and any gc.reflogExpire settings
%
% Not shelling out matters here for two reasons. It is far faster when
% called across ~30 repositories, and MATLAB's system()/! can fail or hang
% outright in some sessions (it does under the MATLAB MCP server, where the
% child process is never even spawned), which would take the whole
% provenance snapshot down with it.
%
% What this CANNOT do without git is determine whether the working tree is
% dirty, or produce a diff - both require hashing the working tree against
% the index. LaBGAScore_prov_snapshot handles that separately and degrades
% gracefully when git is unavailable.
%
%
% *OPTIONS*
%
% * repopath        path to the repository working tree
%
%
% *OUTPUT*
%
% * G               struct with fields
%                   .path, .isrepo, .gitdir
%                   .branch      branch name, or '' when detached
%                   .detached    logical
%                   .commit      full sha of HEAD
%                   .commit_short
%                   .remote      origin url, or ''
%                   .reflog      table: Commit, Time (datetime), Author,
%                                Message - OLDEST FIRST
%                   .reflog_protected  logical, gc.reflogExpire == never
%                   .error       message when something could not be read
%
%
% *SEE ALSO*
%
% LaBGAScore_prov_protect_reflogs, LaBGAScore_prov_snapshot
%
% -------------------------------------------------------------------------
%
% modified by: Lukas Van Oudenhove & Claude Opus 5
%
% date:   KU Leuven, August, 2026
%
% -------------------------------------------------------------------------
%
% LaBGAScore_prov_gitinfo.m         v1.0
%
% last modified: 2026/08/31
%
%
%% INITIALISE
% -------------------------------------------------------------------------

repopath = char(repopath);

G = struct('path', repopath, 'isrepo', false, 'gitdir', '', ...
    'branch', '', 'detached', false, 'commit', '', 'commit_short', '', ...
    'remote', '', 'reflog', [], 'reflog_protected', false, 'error', '');

G.reflog = table(strings(0,1), NaT(0,1), strings(0,1), strings(0,1), ...
    'VariableNames', {'Commit','Time','Author','Message'});


%% LOCATE THE GIT DIRECTORY
% -------------------------------------------------------------------------

gitpath = fullfile(repopath,'.git');

if isfolder(gitpath)
    G.gitdir = gitpath;
elseif isfile(gitpath)
    % worktree or submodule: ".git" is a file containing "gitdir: <path>"
    txt = strtrim(fileread(gitpath));
    tok = regexp(txt,'^gitdir:\s*(.+)$','tokens','once','lineanchors');
    if isempty(tok)
        G.error = 'unreadable .git file';
        return
    end
    gd = strtrim(tok{1});
    if ~startsWith(gd, filesep)
        gd = fullfile(repopath, gd);
    end
    G.gitdir = gd;
else
    G.error = 'not a git repository';
    return
end

G.isrepo = true;


%% RESOLVE HEAD
% -------------------------------------------------------------------------

headfile = fullfile(G.gitdir,'HEAD');

if ~isfile(headfile)
    G.error = 'no HEAD file';
    return
end

head = strtrim(fileread(headfile));
tok = regexp(head,'^ref:\s*(.+)$','tokens','once');

if isempty(tok)
    % detached HEAD: the file holds the sha itself
    G.detached = true;
    G.commit = head;
else
    ref = strtrim(tok{1});
    G.branch = char(extractAfter(string(ref),'refs/heads/'));
    if isempty(G.branch), G.branch = ref; end

    reffile = fullfile(G.gitdir, ref);
    if isfile(reffile)
        G.commit = strtrim(fileread(reffile));
    else
        G.commit = local_packed_ref(G.gitdir, ref);
    end
end

if ~isempty(G.commit) && strlength(G.commit) >= 8
    G.commit_short = G.commit(1:8);
end


%% READ THE REMOTE URL
% -------------------------------------------------------------------------

cfgfile = fullfile(G.gitdir,'config');

if isfile(cfgfile)
    cfg = string(splitlines(string(fileread(cfgfile))));

    insection = false;
    for k = 1:numel(cfg)
        line = strtrim(cfg(k));
        if startsWith(line,'[')
            insection = ~isempty(regexp(line,'^\[remote\s+"origin"\]','once'));
            continue
        end
        if insection
            t = regexp(line,'^url\s*=\s*(.+)$','tokens','once');
            if ~isempty(t), G.remote = strtrim(char(t{1})); break, end
        end
    end

    % gc.reflogExpire - may be written as [gc] reflogExpire = never
    ingc = false;
    for k = 1:numel(cfg)
        line = strtrim(cfg(k));
        if startsWith(line,'[')
            ingc = ~isempty(regexp(line,'^\[gc\]','once'));
            continue
        end
        if ingc
            t = regexp(line,'^reflogExpire\s*=\s*(.+)$','tokens','once','ignorecase');
            if ~isempty(t) && strcmpi(strtrim(char(t{1})),'never')
                G.reflog_protected = true;
            end
        end
    end
end


%% PARSE THE REFLOG
% -------------------------------------------------------------------------
% .git/logs/HEAD format, one entry per line, oldest first:
%   <old-sha> <new-sha> <name> <email> <unixtime> <tz>\t<message>
% The epoch seconds are absolute, so <tz> is informational only.

logfile = fullfile(G.gitdir,'logs','HEAD');

if ~isfile(logfile)
    return
end

lines = splitlines(string(fileread(logfile)));
lines = lines(strlength(strtrim(lines)) > 0);

n = numel(lines);
commit = strings(n,1); time = NaT(n,1); author = strings(n,1); message = strings(n,1);
keep = false(n,1);

for k = 1:n

    parts = split(lines(k), sprintf('\t'));
    meta = parts(1);
    if numel(parts) > 1
        message(k) = strjoin(parts(2:end), sprintf('\t'));
    end

    t = regexp(meta, '^([0-9a-f]{40})\s+([0-9a-f]{40})\s+(.*?)\s+(\d+)\s+([+-]\d{4})\s*$', ...
        'tokens', 'once');

    if isempty(t), continue, end

    commit(k) = t(2);                       % the sha HEAD moved TO
    author(k) = strtrim(t(3));

    % epoch seconds -> local wall-clock time, unzoned, so these compare
    % directly against the file mtimes that dir() returns
    dt = datetime(str2double(t(4)), 'ConvertFrom','posixtime', 'TimeZone','UTC');
    dt.TimeZone = 'local';
    dt.TimeZone = '';
    time(k) = dt;

    keep(k) = true;

end

G.reflog = table(commit(keep), time(keep), author(keep), message(keep), ...
    'VariableNames', {'Commit','Time','Author','Message'});

end


%% LOCAL FUNCTIONS
% -------------------------------------------------------------------------

function sha = local_packed_ref(gitdir, ref)
% look the ref up in .git/packed-refs when it has no loose file

sha = '';

pf = fullfile(gitdir,'packed-refs');
if ~isfile(pf), return, end

lines = splitlines(string(fileread(pf)));

for k = 1:numel(lines)
    line = strtrim(lines(k));
    if startsWith(line,'#') || startsWith(line,'^') || strlength(line) == 0
        continue
    end
    t = regexp(line, '^([0-9a-f]{40})\s+(.+)$', 'tokens', 'once');
    if ~isempty(t) && strcmp(strtrim(char(t(2))), ref)
        sha = char(t(1));
        return
    end
end

end
