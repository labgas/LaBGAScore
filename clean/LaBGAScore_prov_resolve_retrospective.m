function PROVTABLE = LaBGAScore_prov_resolve_retrospective(projdir, varargin)
%
%
% *USAGE*
%
% Reconstructs, for analyses that were ALREADY run, which commit of each
% dependency was checked out at the time. Use it to document existing
% projects, where no provenance was recorded when the scripts ran.
%
% It works because two records survive:
%
% # every html report published by MATLAB carries the run date in a
%   <meta name="DC.date"> tag, the script name in DC.source, the MATLAB
%   version in the generator tag, and the complete script source between
%   its "##### SOURCE BEGIN #####" markers
% # every git clone keeps a reflog in .git/logs/HEAD recording each commit
%   HEAD actually moved to, and when
%
% Intersecting the two gives the commit that was genuinely checked out
% locally when the report was produced. That is a real reconstruction, not
% an inference from commit dates - a commit's own date says nothing about
% when this machine moved to it.
%
% Output is written as SIDECAR files. Existing reports are never modified:
% they are typically already datalad-saved and git-annexed, and rewriting
% them would mean unlocking content just to append a footer.
%
%
% *OPTIONS*
%
% * projdir         project superdataset root, e.g. '/data/proj_cfs'
% * 'subdirs'       which trees to search, default {'secondlevel','firstlevel'}
% * 'repos'         dependency repos to resolve, default: all git repos
%                   under githubrootdir that have a reflog
% * 'githubrootdir' default '/data/master_github_repos'
% * 'depmap'        logical, run a dependency map per report so that dirty
%                   state can be narrowed per script. Default true.
% * 'index'         prebuilt index from LaBGAScore_dep_build_index. Strongly
%                   recommended: without it the index is rebuilt for every
%                   report, which dominates the run time.
% * 'maxdepth'      transitive depth of the per-report dependency map,
%                   default Inf. Lower it to trade precision for speed.
% * 'followrepos'   restrict the transitive walk to these repositories, and
%                   the fallback repository list for artifacts whose
%                   producing script could not be identified.
%                   Recommended: pass the git repositories only. SPM makes
%                   up the bulk of a typical closure but is not version
%                   controlled, so walking it costs time and tells you
%                   nothing about provenance.
% * 'write'         logical, write sidecar files, default true
% * 'print'         logical, default true
%
%
% *OUTPUT*
%
% * PROVTABLE       table, one row per report x dependency, with columns
%                   Model, Report, Script, RunTime, TimeSource, MATLAB,
%                   Repo, Commit, CommitShort, Resolution, SameDayMove,
%                   DirtyStatus, RelevantModifiedFiles, ScriptChanged
%
% Sidecars written per model:
%   <model>/results/notes/provenance_<model>.tsv and .mat
%   <model>/results/html/<report>_provenance.html
% and one project-level provenance_<proj>.tsv at projdir.
%
%
% *RESOLUTION VALUES*
%
% * reflog_exact     the report's timestamp was known to the second and
%                    fell unambiguously between two reflog entries
% * reflog_day       only the run DATE was recoverable (the report is
%                    git-annexed, so its symlink mtime is the annex-add
%                    time, not the run time). End of day is used, which is
%                    conservative. Almost always still unambiguous, since
%                    dependency updates are weeks apart.
% * reflog_day_APPROX  as above, but the repository ALSO moved on that same
%                    calendar day, so the commit given may be one move too
%                    late. SameDayMove is true for these.
% * predates_clone_APPROX  the run is older than this clone's reflog, so
%                    the commit cannot be recovered. Nothing is guessed.
% * no_reflog        the repository has no reflog at all
%
%
% *DIRTY STATUS*
%
% Whether a dependency carried uncommitted edits AT THE TIME cannot be
% recovered - uncommitted work is not timestamped anywhere. What CAN be
% established is whether today's uncommitted edits touch files this script
% reaches, which bounds the uncertainty rather than leaving it open:
%
% * clean_for_this_script  nothing this script reaches is modified today,
%                          so the commit hash is a complete description
%                          unless an edit was made and later reverted
% * dirty_now_relevant     N files this script reaches carry uncommitted
%                          edits today, named in RelevantModifiedFiles
% * unknown                no dependency map available for this report
%
% ScriptChanged is "script_deleted" when the script that produced a report no
% longer exists in the code subdataset. Those runs are still mapped, from the
% source embedded in the report itself.
%
%
% *SEE ALSO*
%
% LaBGAScore_prov_snapshot, LaBGAScore_prov_gitinfo, LaBGAScore_dep_map
%
% -------------------------------------------------------------------------
%
% modified by: Lukas Van Oudenhove & Claude Opus 5
%
% date:   KU Leuven, August, 2026
%
% -------------------------------------------------------------------------
%
% LaBGAScore_prov_resolve_retrospective.m         v1.0
%
% last modified: 2026/08/31
%
%
%% PARSE OPTIONS
% -------------------------------------------------------------------------

p = inputParser;
p.addParameter('subdirs', {'secondlevel','firstlevel'}, @iscell);
p.addParameter('artifacts', {'report','result'}, @iscell);
p.addParameter('repos', {}, @iscell);
p.addParameter('githubrootdir', '/data/master_github_repos', @(x) ischar(x) || isstring(x));
p.addParameter('depmap', true, @islogical);
p.addParameter('index', [], @(x) isempty(x) || isstruct(x));
p.addParameter('maxdepth', Inf, @isnumeric);
p.addParameter('followrepos', {}, @iscell);
p.addParameter('write', true, @islogical);
p.addParameter('print', true, @islogical);
p.parse(varargin{:});
opt = p.Results;

projdir = char(projdir);
githubrootdir = char(opt.githubrootdir);

if ~isfolder(projdir)
    error('projdir %s does not exist', projdir);
end

[~, projname] = fileparts(projdir);


%% BUILD THE REFLOG TIMELINES
% -------------------------------------------------------------------------

if isempty(opt.repos)
    dd = dir(githubrootdir);
    dd = dd([dd.isdir]);
    dd = dd(~startsWith({dd.name},'.'));
    repolist = {dd.name};
else
    repolist = opt.repos;
end

if opt.print
    fprintf('\nReading reflogs from %s\n', githubrootdir);
end

% Repositories to fall back on when an artifact has no dependency map.
% followrepos, when given, is the caller saying which repos matter here.
corerepos = opt.followrepos;
if isempty(corerepos) && ~isempty(opt.repos)
    corerepos = opt.repos;
end

TL = struct('repo',{},'times',{},'commits',{},'gitinfo',{});

for k = 1:numel(repolist)
    G = LaBGAScore_prov_gitinfo(fullfile(githubrootdir, repolist{k}));
    if ~G.isrepo || height(G.reflog) == 0, continue, end
    TL(end+1).repo = repolist{k}; %#ok<AGROW>
    TL(end).times = G.reflog.Time;                  % oldest first
    TL(end).commits = G.reflog.Commit;
    TL(end).gitinfo = G;
end

if opt.print
    fprintf('%d repositories have a usable reflog\n', numel(TL));
end


%% FIND THE PUBLISHED REPORTS
% -------------------------------------------------------------------------

reports = {};

for s = 1:numel(opt.subdirs)
    base = fullfile(projdir, opt.subdirs{s});
    if ~isfolder(base), continue, end
    d = dir(fullfile(base,'**','*.html'));
    d = d(~[d.isdir]);
    for k = 1:numel(d)
        % annexed content lives under .git/annex/objects as a second copy of
        % the same report; including it would invent duplicate runs
        if contains([d(k).folder filesep], [filesep '.git' filesep]), continue, end
        % our own sidecars are html too, and would otherwise be picked up as
        % reports on a second run: <report>_provenance.html next to each
        % report, and provenance_<model>.html for the per-model overview
        if endsWith(d(k).name, '_provenance.html'), continue, end
        if startsWith(d(k).name, 'provenance_'), continue, end
        reports{end+1} = fullfile(d(k).folder, d(k).name); %#ok<AGROW>
    end
end

% Result .mat files matter as much as reports. Only a few scripts per model
% are ever published to html - the prep_ scripts write .mat files instead -
% so reports alone leave most of the pipeline undocumented. A .mat carries an
% embedded "Created on" stamp in its 116-byte header, which is the same kind
% of evidence as a report's DC.date and survives git-annex just as well.

results = {};

if ismember('result', opt.artifacts)
    for s = 1:numel(opt.subdirs)
        base = fullfile(projdir, opt.subdirs{s});
        if ~isfolder(base), continue, end
        d = dir(fullfile(base,'**','results','**','*.mat'));
        d = d(~[d.isdir]);
        for k = 1:numel(d)
            if contains([d(k).folder filesep], [filesep '.git' filesep]), continue, end
            if startsWith(d(k).name, 'provenance_'), continue, end   % our own output
            if strcmpi(d(k).name, 'SPM.mat'), continue, end          % SPM's own design file
            results{end+1} = fullfile(d(k).folder, d(k).name); %#ok<AGROW>
        end
    end
end

if ~ismember('report', opt.artifacts)
    reports = {};
end

artifacts = [reports(:); results(:)];
artifacttype = [repmat("report", numel(reports), 1); repmat("result", numel(results), 1)];

if opt.print
    fprintf('%d published html report(s) and %d result file(s) found under %s\n\n', ...
        numel(reports), numel(results), projdir);
end

if isempty(artifacts)
    PROVTABLE = local_empty_table();
    return
end


%% RESOLVE EACH REPORT
% -------------------------------------------------------------------------

model = strings(0,1); report = strings(0,1); atype = strings(0,1); script = strings(0,1);
mindepth = zeros(0,1);
runtime = strings(0,1); timesource = strings(0,1); matlabver = strings(0,1);
repocol = strings(0,1); commit = strings(0,1); commitshort = strings(0,1);
resolution = strings(0,1); samedaymove = false(0,1);
dirtystatus = strings(0,1); relevantmod = strings(0,1); scriptchanged = strings(0,1);

INFO = struct('report',{},'model',{},'rows',{});
DEPCACHE = containers.Map('KeyType','char','ValueType','any');
SRCCACHE = containers.Map('KeyType','char','ValueType','any');
SCRIPTCACHE = containers.Map('KeyType','char','ValueType','any');

for r = 1:numel(artifacts)

    thisreport = artifacts{r};
    thistype = artifacttype(r);

    if thistype == "report"
        M = local_parse_report(thisreport);
    else
        M = local_parse_resultfile(thisreport);
    end

    thismodel = local_model_of(thisreport, projdir);

    % A result file whose name follows no convention we know can still often
    % be attributed: ask which of this model's own scripts mentions it. That
    % generalises to study-specific outputs (The Decoding Toolbox's res_*.mat,
    % JuSpace's neurotransmitter_*.mat) without guessing at conventions.
    if thistype == "result" && isempty(M.script)
        M.script = local_script_from_sources(projdir, thismodel, thisreport, SRCCACHE);
    end

    % A convention-derived name is the CANlab template's, e.g.
    % prep_3a_run_second_level_regression_and_save. The study's copy is
    % renamed AND usually truncated - cfs_secondlevel_m14_s6_prep_3a_run_regression -
    % so it shares no suffix with the template. Map it onto this model's own
    % script, otherwise a report and the .mat files from the same run end up
    % under two different names and get two provenance pages instead of one.
    if thistype == "result" && ~isempty(M.script)
        mapped = local_study_script(projdir, thismodel, M.script, SCRIPTCACHE);
        if ~isempty(mapped), M.script = mapped; end
    end

    % --- dependency map for this report's script -------------------------

    usedfiles = strings(0,1);
    haddepmap = false;
    mdrepo = strings(0,1); mdval = zeros(0,1);

    if opt.depmap && ~isempty(M.script)
        scriptpath = local_find_script(projdir, M.script);
        % the script may have been deleted or renamed since it ran. The report
        % embeds its full source, so recover it to a temporary file and map
        % that instead - these are exactly the runs whose provenance would
        % otherwise be lost completely
        if isempty(scriptpath) && ~isempty(M.source)
            scriptpath = local_source_to_tempfile(M.script, M.source);
        end
        if ~isempty(scriptpath)
            % several reports can come from the same script, and a map takes
            % tens of seconds, so build each one only once
            if isKey(DEPCACHE, scriptpath)
                cached = DEPCACHE(scriptpath);
                usedfiles = cached.files; mdrepo = cached.mdrepo; mdval = cached.mdval;
                haddepmap = true;
            else
                try
                    DEP = LaBGAScore_dep_map(scriptpath, 'print', false, ...
                        'index', opt.index, 'maxdepth', opt.maxdepth, ...
                        'followrepos', opt.followrepos);
                    usedfiles = DEP.files.File;
                    e = DEP.edges(DEP.edges.Repo ~= "",:);
                    if ~isempty(e)
                        [mdrepo, ~, ix] = unique(e.Repo);
                        mdval = accumarray(ix, e.Depth, [], @min);
                    end
                    DEPCACHE(scriptpath) = struct('files', usedfiles, ...
                        'mdrepo', mdrepo, 'mdval', mdval);
                    haddepmap = true;
                catch
                end
            end
        end
    end

    % --- did the script change since it produced this report? ------------

    changed = "unknown";
    if ~isempty(M.source)
        diskpath = local_find_script(projdir, M.script);
        if isempty(diskpath)
            changed = "script_deleted";
        elseif isfile(diskpath)
            now_src = string(fileread(diskpath));
            if strcmp(local_norm(now_src), local_norm(M.source))
                changed = "no";
            else
                changed = "YES";
            end
        end
    end

    % --- resolve every dependency at that moment -------------------------

    nrows = 0;

    for k = 1:numel(TL)

        [thiscommit, thisres, thissameday] = local_head_at(TL(k), M.runtime, M.timesource);

        % narrow dirty state to files this script actually reaches
        thispath = fullfile(githubrootdir, TL(k).repo);
        inrepo = usedfiles(startsWith(usedfiles, string(thispath) + filesep));

        if ~haddepmap
            % No dependency map for this artifact, so we cannot say which
            % repositories it actually uses. Emitting a row for all ~29 clones
            % would bury the record in repos nothing here touches
            % (BrainNetClass, ExploreASL, the website). Restrict to the ones
            % the caller declared an interest in; the commit resolution is
            % still recorded, only the per-script narrowing is missing.
            if ~isempty(corerepos) && ~ismember(TL(k).repo, corerepos)
                continue
            end
            thisdirty = "unknown"; thisrel = "";
        elseif isempty(inrepo)
            continue                        % this script never touches this repo
        else
            relused = erase(inrepo, string(thispath) + filesep);
            S = LaBGAScore_prov_gitstatus(thispath, 'subset', relused);
            if ~S.ok
                thisdirty = "unknown"; thisrel = "";
            elseif isempty(S.modified) && isempty(S.deleted)
                thisdirty = "clean_for_this_script"; thisrel = "";
            else
                thisdirty = "dirty_now_relevant";
                thisrel = strjoin([S.modified; S.deleted], '; ');
            end
        end

        model(end+1,1) = string(thismodel); %#ok<AGROW>
        report(end+1,1) = string(thisreport); %#ok<AGROW>
        atype(end+1,1) = thistype; %#ok<AGROW>
        md = mdval(mdrepo == string(TL(k).repo));
        if isempty(md), mindepth(end+1,1) = NaN; else, mindepth(end+1,1) = md(1); end %#ok<AGROW>
        script(end+1,1) = string(M.script); %#ok<AGROW>
        runtime(end+1,1) = string(M.runtime, 'yyyy-MM-dd''T''HH:mm:ss'); %#ok<AGROW>
        timesource(end+1,1) = string(M.timesource); %#ok<AGROW>
        matlabver(end+1,1) = string(M.matlab); %#ok<AGROW>
        repocol(end+1,1) = string(TL(k).repo); %#ok<AGROW>
        commit(end+1,1) = string(thiscommit); %#ok<AGROW>
        if strlength(thiscommit) >= 8
            commitshort(end+1,1) = extractBefore(string(thiscommit),9); %#ok<AGROW>
        else
            commitshort(end+1,1) = ""; %#ok<AGROW>
        end
        resolution(end+1,1) = string(thisres); %#ok<AGROW>
        samedaymove(end+1,1) = thissameday; %#ok<AGROW>
        dirtystatus(end+1,1) = thisdirty; %#ok<AGROW>
        relevantmod(end+1,1) = thisrel; %#ok<AGROW>
        scriptchanged(end+1,1) = changed; %#ok<AGROW>

        nrows = nrows + 1;

    end

    INFO(end+1).report = thisreport; %#ok<AGROW>
    INFO(end).model = thismodel;
    INFO(end).rows = nrows;

    if opt.print
        [~, rn, re] = fileparts(thisreport);
        fprintf('  %-56s %s (%s)\n', [rn re], string(M.runtime,'yyyy-MM-dd'), M.timesource);
    end

end

PROVTABLE = table(model, report, atype, script, runtime, timesource, matlabver, ...
    repocol, mindepth, commit, commitshort, resolution, samedaymove, ...
    dirtystatus, relevantmod, scriptchanged, ...
    'VariableNames', {'Model','Artifact','ArtifactType','Script','RunTime','TimeSource','MATLAB', ...
        'Repo','MinDepth','Commit','CommitShort','Resolution','SameDayMove', ...
        'DirtyStatus','RelevantModifiedFiles','ScriptChanged'});

% direct dependencies first, then by how deeply the repo is reached
PROVTABLE = sortrows(PROVTABLE, {'Model','Artifact','MinDepth','Repo'});


%% WRITE SIDECARS
% -------------------------------------------------------------------------

if opt.write && ~isempty(PROVTABLE)
    local_write_sidecars(PROVTABLE, projdir, projname, opt.subdirs, opt.print);
end


%% PRINT SUMMARY
% -------------------------------------------------------------------------

if opt.print
    local_print_summary(PROVTABLE);
end

end


%% LOCAL FUNCTIONS
% -------------------------------------------------------------------------

function M = local_parse_report(f)
% pull run date, script name, MATLAB version and embedded source out of a
% report published by MATLAB

M = struct('script','','runtime',NaT,'timesource','','matlab','','source','');

txt = fileread(f);

t = regexp(txt,'name="DC\.source" content="([^"]+)"','tokens','once');
if ~isempty(t), M.script = char(t{1}); end

t = regexp(txt,'name="generator" content="([^"]+)"','tokens','once');
if ~isempty(t), M.matlab = char(t{1}); end

% the complete script as it was when it ran
t = regexp(txt,'##### SOURCE BEGIN #####(.*?)##### SOURCE END #####','tokens','once');
if ~isempty(t)
    % publish() cannot leave "--" inside an HTML comment, so it substitutes
    % a placeholder; undo that or the source will never match the file
    M.source = strrep(char(t{1}), 'REPLACE_WITH_DASH_DASH', '--');
end

% run date from the DC.date tag; time of day from the file itself, which is
% only meaningful when the two agree - for a git-annexed report the file is
% a symlink whose mtime is the annex-add time, not the run time
dcdate = NaT;
t = regexp(txt,'name="DC\.date" content="([^"]+)"','tokens','once');
if ~isempty(t)
    dcdate = datetime(t{1},'InputFormat','yyyy-MM-dd');
end

d = dir(f);
fmtime = datetime(d.datenum,'ConvertFrom','datenum');

if isnat(dcdate)
    M.runtime = fmtime;
    M.timesource = 'mtime_only';
elseif dateshift(fmtime,'start','day') == dcdate
    M.runtime = fmtime;
    M.timesource = 'DC.date+mtime';
else
    M.runtime = dcdate + hours(23) + minutes(59) + seconds(59);
    M.timesource = 'DC.date+eod';
end

end


function M = local_parse_resultfile(f)
% Run metadata for a result .mat file.
%
% Every MAT-file begins with a 116-byte text header of the form
%   MATLAB 7.3 MAT-file, Platform: GLNXA64, Created on: Fri Aug 28 21:40:23 2026
% That stamp is written into the file, so unlike the filesystem mtime it
% survives git-annex, copying and datalad. It is the .mat equivalent of the
% DC.date tag in a published report, and carries a time of day as well.
%
% Caveat for files built with save(...,'-append'): the header keeps the
% ORIGINAL creation time while the mtime moves to the last append. Both are
% recorded, and the later of the two is used, since that is when the file
% last had content written to it.

M = struct('script','','runtime',NaT,'timesource','','matlab','','source','');

created = NaT;

fid = fopen(f,'r');
if fid > 0
    hdr = fread(fid, 116, '*char')';
    fclose(fid);
    % match the date pattern exactly: the 116-byte header also carries
    % trailing text such as "HDF5 schema 1.00", which a lazy .+? would
    % swallow. The day may be space-padded ("Aug  5"), hence \s+.
    tok = regexp(hdr, ...
        'Created on:\s*(\w{3}\s+\w{3}\s+\d{1,2}\s+\d{2}:\d{2}:\d{2}\s+\d{4})', ...
        'tokens', 'once');
    if ~isempty(tok)
        try
            created = datetime(regexprep(strtrim(tok{1}), '\s+', ' '), ...
                'InputFormat','eee MMM d HH:mm:ss yyyy', 'Locale','en_US');
        catch
        end
    end
    v = regexp(hdr, 'MATLAB (\d+\.\d+) MAT-file', 'tokens', 'once');
    if ~isempty(v), M.matlab = ['MAT-file v' v{1}]; end
end

d = dir(f);
mtime = datetime(d.datenum, 'ConvertFrom','datenum');

% An annexed file is a symlink into .git/annex/objects whose mtime is the
% annex-add time, not the write time; the embedded stamp is all we can trust.
isannexed = false;
try
    isannexed = ~isempty(regexp(char(java.io.File(f).getCanonicalPath()), ...
        [filesep '\.git' filesep 'annex' filesep], 'once'));
catch
end

if isnat(created)
    M.runtime = mtime;
    M.timesource = 'mtime_only';
elseif isannexed
    M.runtime = created;
    M.timesource = 'mat_header';
elseif mtime > created
    M.runtime = mtime;                     % appended after creation
    M.timesource = 'mat_header+append';
else
    M.runtime = created;
    M.timesource = 'mat_header';
end

M.script = local_script_for_mat(f);

end


function name = local_study_script(projdir, modelname, canonical, cache)
% Find this model's copy of a CANlab template script.
%
% Copies are renamed <proj>_<level>_m<M>_s<N>_<something>, and the trailing
% part is usually shortened, so neither an exact nor a suffix match works.
% Compare the token sets instead, requiring the step designator (prep_3a, c2a,
% a2 ...) to agree so that prep_3 and prep_3a can never be confused.

name = '';

if isKey(cache, modelname)
    S = cache(modelname);
else
    d = dir(fullfile(projdir, 'code', '**', modelname, '*.m'));
    d = d(~[d.isdir]);
    S = {d.name};
    cache(modelname) = S;
end

if isempty(S), return, end

[~, canon] = fileparts(canonical);
ctok = local_tokens(canon);
ckey = local_stepkey(canon);
if isempty(ctok), return, end

best = ''; bestscore = 0;

for k = 1:numel(S)
    [~, b] = fileparts(S{k});
    % strip the <proj>_<level>_m<M>_s<N>_ prefix to get the comparable part
    rem = regexp(b, '^[A-Za-z0-9\-]+_(?:secondlevel|firstlevel|prep|pet|mrs)_(?:m\d+_)?(?:s\d+[a-z]?_)?(.+)$', ...
        'tokens', 'once');
    if isempty(rem), continue, end
    rtok = local_tokens(rem{1});
    if ~strcmp(local_stepkey(rem{1}), ckey), continue, end
    score = numel(intersect(rtok, ctok)) / numel(ctok);
    if score > bestscore, best = b; bestscore = score; end
end

if bestscore >= 0.4
    name = [best '.m'];
end

end


function t = local_tokens(str)
t = split(string(lower(str)), '_');
t = cellstr(t(strlength(t) > 0));
end


function k = local_stepkey(str)
% leading step designator: prep_3a, prep_1b, c2a, a2, h1, d10 ...
m = regexp(lower(str), '^(prep_?\d*[a-z]?|[a-z]\d*[a-z]?)(?:_|$)', 'tokens', 'once');
if isempty(m), k = ''; else, k = strrep(m{1}, '_', ''); end
end


function name = local_script_from_sources(projdir, modelname, matfile, cache)
% Attribute a result file to a script by looking for its name in the model's
% own sources. Deliberately conservative: only when EXACTLY ONE script
% mentions it, so a shared filename never produces a confident wrong answer.
%
% Result filenames are often assembled at run time, e.g.
%   ['regression_stats_and_maps_' mygroupnamefield '_' scaling_string '.mat']
% so the full basename may appear nowhere in the source. Progressively
% shorter underscore-delimited prefixes are tried, longest first, which
% recovers the static part of such a name.

name = '';

[~, base] = fileparts(matfile);

% scripts of this model, read once and cached
if isKey(cache, modelname)
    S = cache(modelname);
else
    d = dir(fullfile(projdir, 'code', '**', modelname, '*.m'));
    d = d(~[d.isdir]);
    S = struct('name', {}, 'text', {});
    for k = 1:numel(d)
        try
            S(end+1).name = d(k).name; %#ok<AGROW>
            S(end).text = lower(fileread(fullfile(d(k).folder, d(k).name)));
        catch
        end
    end
    cache(modelname) = S;
end

if isempty(S), return, end

parts = strsplit(lower(base), '_');

for n = numel(parts):-1:2
    needle = strjoin(parts(1:n), '_');
    if strlength(needle) < 6, break, end
    hit = find(contains({S.text}, needle));
    if numel(hit) == 1
        [~, b] = fileparts(S(hit).name);
        name = [b '.m'];        % keep the extension, as DC.source and
        return                  % local_study_script do, or the same script
                                % appears under two names and splits its group
    elseif numel(hit) > 1
        return                              % ambiguous: say nothing
    end
end

end


function name = local_script_for_mat(f)
% Best-effort guess at which script wrote a result file, from the CANlab
% naming convention. Returns '' when there is no reliable mapping - the
% commit resolution below does not depend on knowing the script, it only
% gets more specific when the script IS known, because then the call graph
% can narrow which uncommitted changes matter.

[~, base] = fileparts(f);
base = lower(base);

if startsWith(base, 'image_names_and_setup')
    name = 'prep_1_set_conditions_contrasts_colors';
elseif startsWith(base, 'data_objects')
    name = 'prep_2_load_image_data_and_save';
elseif startsWith(base, 'contrast_data_objects')
    name = 'prep_3_calc_univariate_contrast_maps_and_save';
elseif startsWith(base, 'regression_stats_and_maps') || ...
       startsWith(base, 'parcelwise_stats_and_maps') || ...
       startsWith(base, 'mvpa_stats_and_maps') || ...
       startsWith(base, 'roi_stats')
    name = 'prep_3a_run_second_level_regression_and_save';
elseif startsWith(base, 'svm_stats')
    name = 'prep_3c_run_SVMs_on_contrasts_masked';
else
    name = '';                              % study-specific or unknown
end

end


function [commit, resolution, sameday] = local_head_at(tl, t, timesource)
% the commit HEAD had moved to as of time t, from this repo's reflog

commit = ''; sameday = false;

idx = find(tl.times <= t, 1, 'last');

if isempty(idx)
    resolution = 'predates_clone_APPROX';
    return
end

commit = char(tl.commits(idx));

% did this repository also move on the same calendar day? if so, and we
% only know the run date, the commit above could be one move too late
sameday = any(dateshift(tl.times,'start','day') == dateshift(t,'start','day'));

% these sources carry a time of day, so the resolution is to the second;
% only a report dated from DC.date alone is limited to day precision
exactsources = {'DC.date+mtime','mat_header','mat_header+append','mtime_only'};

if ismember(timesource, exactsources)
    resolution = 'reflog_exact';
elseif sameday
    resolution = 'reflog_day_APPROX';
else
    resolution = 'reflog_day';
end

end


function m = local_model_of(reportpath, projdir)
% .../secondlevel/model_14_IOM_pos_neg/results/html/x.html -> model_14_...

rel = erase(reportpath, [projdir filesep]);
parts = strsplit(rel, filesep);
if numel(parts) >= 2, m = parts{2}; else, m = 'unknown'; end

end


function fp = local_find_script(projdir, scriptname)
% locate the script that produced a report, in the project code subdataset

fp = '';
if isempty(scriptname), return, end

[~, base] = fileparts(scriptname);

d = dir(fullfile(projdir,'code','**',[base '.m']));
d = d(~[d.isdir]);
if ~isempty(d)
    fp = fullfile(d(1).folder, d(1).name);
    return
end

% Names taken from the CANlab convention (prep_2_load_image_data_and_save)
% never match a study copy directly, because copies are renamed
% <proj>_secondlevel_m<M>_s<N>_<template>. Match on that suffix so the
% dependency map analyses the study's own adapted script rather than the
% generic template it came from.
d = dir(fullfile(projdir,'code','**',['*_' base '.m']));
d = d(~[d.isdir]);
if ~isempty(d)
    fp = fullfile(d(1).folder, d(1).name);
    return
end

w = which(base);
if ~isempty(w) && isfile(w), fp = w; end

end


function fp = local_source_to_tempfile(scriptname, source)
% write a report's embedded source out so it can be parsed like any file

fp = '';
[~, base] = fileparts(char(scriptname));
if isempty(base), return, end

fp = fullfile(tempdir, ['labgascore_prov_' base '.m']);

fid = fopen(fp, 'w');
if fid < 0, fp = ''; return, end
fwrite(fid, source);
fclose(fid);

end


function s = local_norm(txt)
% compare sources ignoring line-ending and trailing-whitespace churn
s = regexprep(char(txt), '\r\n?', '\n');
s = regexprep(s, '[ \t]+\n', '\n');
s = strtrim(s);
end


function local_write_sidecars(T, projdir, projname, subdirs, doprint)

models = unique(T.Model);

for k = 1:numel(models)

    sub = T(T.Model == models(k),:);
    % derive the model directory by name: result files and html reports sit at
    % different depths, so walking up a fixed number of levels is not safe
    first = char(sub.Artifact(1));
    ix = strfind(first, [filesep char(models(k)) filesep]);
    modeldir = first(1 : ix(1) + numel(char(models(k))));

    % char throughout: fullfile returns a string when given one, and
    % [string 'suffix'] builds a two-element array rather than concatenating
    % Second-level models keep everything under results/; first-level models
    % have no such directory - subjects sit directly under the model with
    % their reports in sub-*/diagnostics/. Follow whichever layout is already
    % there rather than inventing a results/ tree where the workflow has none.
    if isfolder(fullfile(modeldir,'results'))
        notesdir = char(fullfile(modeldir,'results','notes'));
        htmldir  = char(fullfile(modeldir,'results','html'));
    else
        notesdir = char(fullfile(modeldir,'provenance'));
        htmldir  = notesdir;
    end
    if ~isfolder(notesdir), mkdir(notesdir); end

    stem = char(fullfile(notesdir, sprintf('provenance_%s', models(k))));
    writetable(sub, [stem '.tsv'], 'FileType','text', 'Delimiter','\t');
    PROV = sub; save([stem '.mat'], 'PROV', '-v7.3');

    % one small standalone page per report, next to but never touching it
    if ~isfolder(htmldir), mkdir(htmldir); end

    % ONE page per script run, not per artifact. A single script typically
    % produces both an html report and one or more .mat files; giving each its
    % own page says the same thing several times. Artifacts are grouped by the
    % script that produced them, and the page lists them together, strongest
    % evidence first. Artifacts whose script could not be identified cannot be
    % grouped, so they keep a page of their own.
    grp = strings(height(sub),1);
    for i = 1:height(sub)
        % The group is one RUN of one script, not one script. At second level
        % a script runs once per model, so script name alone is the right key
        % and it merges the report with the .mat files from the same run. At
        % first level the same script runs once per SUBJECT - 135 of them in
        % proj_cfs, spanning five different CanlabCore commits - so the
        % subject has to be part of the key or every subject collapses into a
        % single page and the per-subject versions are lost.
        subj = regexp(char(sub.Artifact(i)), [filesep '(sub-[^' filesep ']+)' filesep], ...
            'tokens', 'once');
        if isempty(subj), subjkey = ""; else, subjkey = "|" + string(subj{1}); end

        if strlength(sub.Script(i)) > 0
            grp(i) = "script:" + sub.Script(i) + subjkey;
        else
            grp(i) = "artifact:" + sub.Artifact(i);
        end
    end

    ugrp = unique(grp);

    for i = 1:numel(ugrp)

        rows = sub(grp == ugrp(i),:);

        % name the page after the report when there is one, so it sits next to
        % the report and is found by anyone reading it; otherwise after the
        % script, in the html directory
        isrep = rows.ArtifactType == "report";
        if any(isrep)
            repfile = rows.Artifact(find(isrep,1));
            [rdir, rname] = fileparts(repfile);
            outfile = char(fullfile(char(rdir), [char(rname) '_provenance.html']));
        elseif strlength(rows.Script(1)) > 0
            [~, sname] = fileparts(rows.Script(1));
            outfile = char(fullfile(htmldir, [char(sname) '_provenance.html']));
        else
            [~, rname, rext] = fileparts(rows.Artifact(1));
            outfile = char(fullfile(htmldir, ...
                [char(rname) char(strrep(rext,'.','_')) '_provenance.html']));
        end

        local_write_html(outfile, rows);

    end

    % and one overview page per model linking them together
    local_write_model_html(char(fullfile(htmldir, ...
        sprintf('provenance_%s.html', models(k)))), sub, char(models(k)));

end

% Summary at the root of each tree scanned (secondlevel/, firstlevel/)
% rather than at the superdataset root: it describes that tree, and belongs
% in the subdataset whose contents it documents.
written = strings(0,1);

for s = 1:numel(subdirs)
    base = fullfile(projdir, subdirs{s});
    if ~isfolder(base), continue, end
    sub = T(startsWith(T.Artifact, string(base) + filesep), :);
    if isempty(sub), continue, end
    out = char(fullfile(base, sprintf('provenance_%s_%s.tsv', projname, subdirs{s})));
    writetable(sub, out, 'FileType','text', 'Delimiter','\t');
    written(end+1,1) = string(out); %#ok<AGROW>
end

if doprint
    fprintf('\nSidecars written for %d model(s). Summaries:\n', numel(models));
    for k = 1:numel(written), fprintf('  %s\n', written(k)); end
end

end


function local_write_html(outfile, rows)
% One page per script run, covering every artifact that run produced.
%
% EVIDENCE HIERARCHY. A script run usually leaves several dated traces, and
% they are not equally good. Strongest first:
%
%   mat_header         "Created on:" inside the .mat file's 116-byte header.
%                      Embedded, carries a time of day, and survives
%                      git-annex, copying and datalad untouched.
%   mat_header+append  same stamp, but the file was appended to later, so the
%                      timestamp used is the later write
%   DC.date+mtime      the report's own run DATE, with time of day taken from
%                      the file's mtime because the two agree
%   DC.date+eod        only the run date was recoverable - the report is
%                      git-annexed, so its mtime is the annex-add time. End of
%                      day is used, which is conservative
%   mtime_only         no embedded date at all; the filesystem is all there is
%
% The best-evidenced artifact drives the dependency table. The others are
% listed with their own resolved commits, so if they disagree - which happens
% when a file is appended to by a later script - that is visible rather than
% averaged away.

nl = newline;

order = ["mat_header","mat_header+append","DC.date+mtime","DC.date+eod","mtime_only"];

arts = unique(rows.Artifact, 'stable');
rank = zeros(numel(arts),1);
for k = 1:numel(arts)
    ts = rows.TimeSource(find(rows.Artifact == arts(k), 1));
    ix = find(order == ts, 1);
    if isempty(ix), rank(k) = numel(order) + 1; else, rank(k) = ix; end
end
[~, ord] = sort(rank);
arts = arts(ord);

primary = rows(rows.Artifact == arts(1),:);

s = ['<!DOCTYPE html><html><head><meta charset="utf-8">' nl];
s = [s '<meta name="viewport" content="width=device-width, initial-scale=1">' nl];
s = [s '<title>Provenance: ' char(primary.Script(1)) '</title>' nl];
s = [s '<style>body{font-family:sans-serif;margin:2em;max-width:64em}' ...
     'table{border-collapse:collapse;font-size:90%;width:100%}' ...
     'th,td{padding:4px 8px;border-bottom:1px solid #ddd;text-align:left;vertical-align:top}' ...
     'th{border-bottom:1px solid #999}code{font-size:95%}' ...
     '.warn{color:#b00}.muted{color:#666}.scroll{overflow-x:auto}' ...
     '</style></head><body>' nl];

s = [s '<h1>Provenance (reconstructed)</h1>' nl];

if strlength(primary.Script(1)) > 0
    s = [s '<p><b><code>' char(primary.Script(1)) '</code></b>'];
    % Only call the attribution inferred when nothing in the group is a
    % published report. A report names its script in DC.source, so if one is
    % present here the identity is recorded, not guessed - even when a .mat
    % file carries the better timestamp and is therefore the primary artifact.
    if ~any(rows.ArtifactType == "report")
        s = [s ' <span class="muted">(script inferred from the file name)</span>'];
    end
    if primary.ScriptChanged(1) == "YES"
        s = [s '<br><span class="warn">This script has changed since it produced these ' ...
             'files: the source embedded in the report differs from the file now in ' ...
             'the code subdataset.</span>'];
    elseif primary.ScriptChanged(1) == "script_deleted"
        s = [s '<br><span class="warn">This script no longer exists in the code ' ...
             'subdataset. It was recovered from the source embedded in the report.</span>'];
    end
    s = [s '</p>' nl];
else
    [~, an, ae] = fileparts(arts(1));
    s = [s '<p><b>' char(an) char(ae) '</b><br>' ...
         '<span class="muted">producing script could not be identified</span></p>' nl];
end

% --- what this run produced -------------------------------------------------

s = [s '<h2>Files</h2>' nl];
s = [s '<div class="scroll"><table>' nl];
s = [s '<tr><th>File</th><th>Type</th><th>Written</th><th>Evidence</th>' ...
     '<th>CanlabCore</th></tr>' nl];

ccommits = strings(0,1);

for k = 1:numel(arts)
    rk = rows(rows.Artifact == arts(k),:);
    [~, an, ae] = fileparts(arts(k));
    cc = rk(rk.Repo == "CanlabCore",:);
    if isempty(cc)
        cstr = '<span class="muted">n/a</span>';
    else
        cstr = ['<code>' char(cc.CommitShort(1)) '</code>'];
        ccommits(end+1,1) = cc.CommitShort(1); %#ok<AGROW>
    end
    if k == 1
        marker = ' <span class="muted">&larr; strongest evidence</span>';
    else
        marker = '';
    end
    s = [s sprintf('<tr><td>%s%s</td><td>%s</td><td>%s</td><td><code>%s</code></td><td>%s</td></tr>%s', ...
        [char(an) char(ae)], marker, char(rk.ArtifactType(1)), char(rk.RunTime(1)), ...
        char(rk.TimeSource(1)), cstr, nl)]; %#ok<AGROW>
end

s = [s '</table></div>' nl];

if numel(unique(ccommits)) > 1
    s = [s '<p class="warn">These files were not all written against the same ' ...
         'CanlabCore commit. That happens when a file is appended to by a later ' ...
         'script; the dependency table below describes the strongest-evidenced ' ...
         'file only.</p>' nl];
end

% --- dependencies, from the best-evidenced artifact -------------------------

s = [s '<h2>Dependencies</h2>' nl];
s = [s '<p class="muted" style="font-size:90%">Resolved by intersecting the time above ' ...
     'with the git reflog of each dependency. The files themselves were not modified.</p>' nl];

s = [s '<div class="scroll"><table>' nl];
s = [s '<tr><th>Repository</th><th>Called at</th><th>Commit</th>' ...
     '<th>Resolution</th><th>Uncommitted changes</th></tr>' nl];

for k = 1:height(primary)
    if primary.DirtyStatus(k) == "dirty_now_relevant"
        dirty = ['<span class="warn">' char(strrep(primary.RelevantModifiedFiles(k),'; ','<br>')) '</span>'];
    elseif primary.DirtyStatus(k) == "clean_for_this_script"
        dirty = '<span class="muted">none affecting this script</span>';
    else
        dirty = '<span class="muted">unknown</span>';
    end
    if isnan(primary.MinDepth(k))
        depth = '<span class="muted">-</span>';
    elseif primary.MinDepth(k) == 1
        depth = 'directly';
    else
        depth = sprintf('depth %d', primary.MinDepth(k));
    end
    s = [s sprintf(['<tr><td>%s</td><td>%s</td><td><code>%s</code></td>' ...
                    '<td>%s</td><td>%s</td></tr>%s'], ...
        char(primary.Repo(k)), depth, char(primary.CommitShort(k)), ...
        char(primary.Resolution(k)), dirty, nl)]; %#ok<AGROW>
end

s = [s '</table></div>' nl];

s = [s '<p class="muted" style="font-size:85%">Whether a dependency carried ' ...
     'uncommitted edits <i>at the time</i> cannot be recovered - uncommitted work is ' ...
     'not timestamped. "none affecting this script" means nothing this script reaches ' ...
     'is modified <i>today</i>. "Called at" is the shortest path from the script to ' ...
     'that repository: depth 1 is a direct call, deeper entries are reached through ' ...
     'another dependency. Full detail for every artifact in this model is in ' ...
     '<code>results/notes/</code>.</p>' nl];
s = [s '</body></html>' nl];

fid = fopen(outfile,'w');
if fid < 0, return, end
fwrite(fid, s);
fclose(fid);

end


function local_write_model_html(outfile, T, modelname)
% One page per model covering every artifact it contains - published reports
% and result files alike. Only a few scripts per model are ever published,
% so a per-report page leaves most of the pipeline invisible; this is the
% view that shows the whole model.

nl = newline;

s = ['<!DOCTYPE html><html><head><meta charset="utf-8">' nl];
s = [s '<meta name="viewport" content="width=device-width, initial-scale=1">' nl];
s = [s '<title>Provenance: ' modelname '</title>' nl];
s = [s '<style>body{font-family:sans-serif;margin:2em;max-width:70em}' ...
     'table{border-collapse:collapse;font-size:90%;width:100%}' ...
     'th,td{padding:4px 8px;border-bottom:1px solid #ddd;text-align:left;vertical-align:top}' ...
     'th{border-bottom:1px solid #999}code{font-size:95%}' ...
     '.warn{color:#b00}.muted{color:#666}.scroll{overflow-x:auto}' ...
     '</style></head><body>' nl];

s = [s '<h1>Provenance: ' modelname '</h1>' nl];
s = [s '<p class="muted">Reconstructed by LaBGAScore_prov_resolve_retrospective. ' ...
     'Every artifact in this model is listed, not only the scripts that were ' ...
     'published to html: result <code>.mat</code> files carry an embedded ' ...
     'creation stamp, which dates them just as a report''s DC.date does. ' ...
     'Existing files were not modified.</p>' nl];

arts = unique(T.Artifact, 'stable');

s = [s '<div class="scroll"><table>' nl];
s = [s '<tr><th>Artifact</th><th>Type</th><th>Written</th><th>Script</th>' ...
     '<th>CanlabCore</th><th>Uncommitted changes</th></tr>' nl];

for k = 1:numel(arts)

    rows = T(T.Artifact == arts(k),:);
    [~, nm, ext] = fileparts(arts(k));

    cc = rows(rows.Repo == "CanlabCore",:);
    if isempty(cc)
        commit = '<span class="muted">n/a</span>'; dirty = '';
    else
        commit = ['<code>' char(cc.CommitShort(1)) '</code>'];
        if cc.DirtyStatus(1) == "dirty_now_relevant"
            dirty = ['<span class="warn">' char(strrep(cc.RelevantModifiedFiles(1),'; ','<br>')) '</span>'];
        elseif cc.DirtyStatus(1) == "clean_for_this_script"
            dirty = '<span class="muted">none affecting this script</span>';
        else
            dirty = '<span class="muted">unknown</span>';
        end
    end

    scriptname = char(rows.Script(1));
    if isempty(scriptname)
        scriptname = '<span class="muted">not identifiable</span>';
    elseif rows.ArtifactType(1) == "result"
        scriptname = [scriptname ' <span class="muted">(inferred)</span>'];
    end
    if rows.ScriptChanged(1) == "YES"
        scriptname = [scriptname ' <span class="warn">(changed since)</span>'];
    elseif rows.ScriptChanged(1) == "script_deleted"
        scriptname = [scriptname ' <span class="warn">(deleted since)</span>'];
    end

    s = [s sprintf(['<tr><td>%s</td><td>%s</td><td>%s<br><span class="muted">%s</span></td>' ...
                    '<td>%s</td><td>%s</td><td>%s</td></tr>%s'], ...
        [char(nm) char(ext)], char(rows.ArtifactType(1)), ...
        char(rows.RunTime(1)), char(rows.TimeSource(1)), ...
        scriptname, commit, dirty, nl)]; %#ok<AGROW>

end

s = [s '</table></div>' nl];

s = [s '<p class="muted" style="font-size:85%">Whether a dependency carried ' ...
     'uncommitted edits <i>at the time</i> cannot be recovered - uncommitted work ' ...
     'is not timestamped. "none affecting this script" means nothing the script ' ...
     'reaches is modified <i>today</i>. Script names for result files are inferred ' ...
     'from the CANlab naming convention; the commit resolution does not depend on ' ...
     'them. Full detail for every dependency, not just CanlabCore, is in ' ...
     '<code>results/notes/provenance_' modelname '.tsv</code>.</p>' nl];
s = [s '</body></html>' nl];

fid = fopen(outfile,'w');
if fid < 0, return, end
fwrite(fid, s);
fclose(fid);

end


function local_print_summary(T)

fprintf('\n%d artifact x dependency rows (%d reports, %d result files)\n', height(T), ...
    numel(unique(T.Artifact(T.ArtifactType == "report"))), ...
    numel(unique(T.Artifact(T.ArtifactType == "result"))));

fprintf('\nResolution:\n');
[c, v] = groupcounts(T.Resolution);
for k = 1:numel(v), fprintf('  %-24s %d\n', v(k), c(k)); end

fprintf('\nDirty status:\n');
[c, v] = groupcounts(T.DirtyStatus);
for k = 1:numel(v), fprintf('  %-24s %d\n', v(k), c(k)); end

flagged = T(T.DirtyStatus == "dirty_now_relevant",:);
if ~isempty(flagged)
    fprintf('\n*** %d row(s) where uncommitted edits touch files the script uses ***\n', height(flagged));
    u = unique(flagged(:,{'Repo','RelevantModifiedFiles'}));
    for k = 1:height(u)
        fprintf('  %s: %s\n', u.Repo(k), u.RelevantModifiedFiles(k));
    end
end

changed = unique(T.Artifact(T.ScriptChanged == "YES"));
if ~isempty(changed)
    fprintf('\n%d report(s) whose script has since changed:\n', numel(changed));
    for k = 1:numel(changed)
        [~,n,e] = fileparts(changed(k));
        fprintf('  %s\n', [char(n) char(e)]);
    end
end

fprintf('\n');

end


function T = local_empty_table()
e = strings(0,1);
T = table(e,e,e,e,e,e,e,e,zeros(0,1),e,e,e,false(0,1),e,e,e, ...
    'VariableNames', {'Model','Artifact','ArtifactType','Script','RunTime','TimeSource','MATLAB', ...
        'Repo','MinDepth','Commit','CommitShort','Resolution','SameDayMove', ...
        'DirtyStatus','RelevantModifiedFiles','ScriptChanged'});
end
