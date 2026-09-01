function R = LaBGAScore_dep_report(repodir, varargin)
%
%
% *USAGE*
%
% Generates the dependency documentation for a repository: which external
% functions each script calls, which repository each of those lives in, and
% where static analysis could not give a definite answer.
%
% Writes three files into the repository root, all GENERATED - never edit
% them by hand, re-run this instead:
%
%   dependencies.tsv    one row per call edge, the machine-readable form
%   dependencies.yml    the same keyed by script, for the LaBGAS website
%   DEPENDENCIES.md     the rendered document, readable on GitHub
%
% Direct calls only by default. That is the useful granularity for
% documentation: the transitive closure of a second-level script runs to
% ~2800 files and says little more than "it uses CanlabCore and SPM". The
% full closure is what LaBGAScore_prov_snapshot uses, for a different
% purpose - deciding whether an uncommitted edit can have affected a run.
%
%
% *OPTIONS*
%
% * repodir         repository to document
% * 'index'         prebuilt index from LaBGAScore_dep_build_index
% * 'maxdepth'      default 1, i.e. direct calls
% * 'files'         explicit cellstr of files to document, instead of every
%                   .m file under repodir. Used for CANlab_help_examples,
%                   where only a subset of the templates is actively used.
% * 'reponame'      name for the document title, default: folder name
% * 'subtitle'      one line under the title
% * 'outdir'        where to write, default repodir
% * 'print'         default true
%
%
% *OUTPUT*
%
% * R               struct with .edges, .scripts, .files_written
%
%
% *SEE ALSO*
%
% LaBGAScore_dep_map, LaBGAScore_dep_build_index, LaBGAScore_prov_snapshot
%
% -------------------------------------------------------------------------
%
% modified by: Lukas Van Oudenhove & Claude Opus 5
%
% date:   KU Leuven, August, 2026
%
% -------------------------------------------------------------------------
%
% LaBGAScore_dep_report.m         v1.0
%
% last modified: 2026/08/31
%
%
%% PARSE OPTIONS
% -------------------------------------------------------------------------

p = inputParser;
p.addParameter('index', [], @(x) isempty(x) || isstruct(x));
p.addParameter('maxdepth', 1, @isnumeric);
p.addParameter('files', {}, @iscell);
p.addParameter('reponame', '', @(x) ischar(x) || isstring(x));
p.addParameter('subtitle', '', @(x) ischar(x) || isstring(x));
p.addParameter('regenhint', '', @(x) ischar(x) || isstring(x));
p.addParameter('outdir', '', @(x) ischar(x) || isstring(x));
p.addParameter('print', true, @islogical);
p.parse(varargin{:});
opt = p.Results;

repodir = char(repodir);
if isempty(opt.outdir), outdir = repodir; else, outdir = char(opt.outdir); end
if isempty(opt.reponame), [~, reponame] = fileparts(repodir); else, reponame = char(opt.reponame); end

if isempty(opt.index)
    IDX = LaBGAScore_dep_build_index('print', opt.print);
else
    IDX = opt.index;
end

if isempty(opt.files)
    d = dir(fullfile(repodir,'**','*.m'));
    d = d(~[d.isdir]);
    files = arrayfun(@(x) fullfile(x.folder,x.name), d, 'UniformOutput', false);
else
    files = opt.files;
end

files = files(:);


%% MAP EVERY FILE
% -------------------------------------------------------------------------

if opt.print
    fprintf('\nMapping %d files in %s\n', numel(files), reponame);
end

alledges = table();

for k = 1:numel(files)
    try
        D = LaBGAScore_dep_map(files{k}, 'index', IDX, 'maxdepth', opt.maxdepth, 'print', false);
        alledges = [alledges; D.edges]; %#ok<AGROW>
    catch ME
        warning('LaBGAScore_dep_report:map','%s: %s', files{k}, ME.message);
    end
    if opt.print && mod(k,25) == 0
        fprintf('  %d/%d\n', k, numel(files));
    end
end

% only external dependencies are interesting here: drop calls that resolve
% back into the repository being documented, and MathWorks' own functions
ext = alledges(alledges.CalleeFile ~= "" & ...
               alledges.Confidence ~= "toolbox" & ...
               ~startsWith(alledges.CalleeFile, string(repodir)), :);

R.edges = ext;
R.allEdges = alledges;


%% PER-SCRIPT SUMMARY
% -------------------------------------------------------------------------

[uf, ~, ~] = unique(string(files));

sname = strings(0,1); sdomain = strings(0,1); srepos = strings(0,1);
sndirect = zeros(0,1); scaveats = zeros(0,1); srel = strings(0,1);

for k = 1:numel(uf)

    [~, base] = fileparts(uf(k));
    rows = ext(ext.CallerFile == uf(k),:);
    allrows = alledges(alledges.CallerFile == uf(k),:);

    rel = erase(uf(k), string(repodir) + filesep);
    parts = split(rel, filesep);
    if numel(parts) > 1, dom = parts(1); else, dom = "(root)"; end

    if isempty(rows)
        repos = "";
    else
        repos = strjoin(unique(rows.Repo)', ', ');
    end

    sname(end+1,1) = base; %#ok<AGROW>
    sdomain(end+1,1) = dom; %#ok<AGROW>
    srepos(end+1,1) = repos; %#ok<AGROW>
    sndirect(end+1,1) = numel(unique(rows.Callee)); %#ok<AGROW>
    scaveats(end+1,1) = sum(ismember(allrows.Confidence, ...
        ["ambiguous","dotcall","unparseable","dynamic"])); %#ok<AGROW>
    srel(end+1,1) = rel; %#ok<AGROW>

end

R.scripts = table(sname, srel, sdomain, srepos, sndirect, scaveats, ...
    'VariableNames', {'Script','Path','Domain','DependsOn','NDirect','NCaveats'});
R.scripts = sortrows(R.scripts, {'Domain','Script'});


%% WRITE THE ARTEFACTS
% -------------------------------------------------------------------------

R.files_written = strings(0,1);

if ~isfolder(outdir), mkdir(outdir); end

tsv = fullfile(outdir,'dependencies.tsv');
writetable(ext, tsv, 'FileType','text', 'Delimiter','\t');
R.files_written(end+1,1) = string(tsv);

yml = fullfile(outdir,'dependencies.yml');
local_write_yaml(yml, R, reponame, IDX);
R.files_written(end+1,1) = string(yml);

md = fullfile(outdir,'DEPENDENCIES.md');
regenhint = char(opt.regenhint);
if isempty(regenhint)
    if isempty(opt.files)
        regenhint = sprintf('LaBGAScore_dep_report(''%s'')', repodir);
    else
        % a file list was supplied, so the bare call would document something
        % else entirely - say so rather than print a command that misleads
        regenhint = sprintf(['LaBGAScore_dep_report(''%s'', ''files'', <the %d files ' ...
            'documented here>, ''outdir'', ''%s'')'], repodir, numel(files), outdir);
    end
end
local_write_md(md, R, reponame, char(opt.subtitle), regenhint, opt.maxdepth);
R.files_written(end+1,1) = string(md);

if opt.print
    fprintf('\n%d files, %d external call edges, %d distinct dependencies\n', ...
        numel(files), height(ext), numel(unique(ext.Callee)));
    fprintf('Written:\n');
    for k = 1:numel(R.files_written), fprintf('  %s\n', R.files_written(k)); end
    fprintf('\n');
end

end


%% LOCAL FUNCTIONS
% -------------------------------------------------------------------------

function local_write_md(mdfile, R, reponame, subtitle, regenhint, maxdepth)

nl = newline;
fid = fopen(mdfile,'w');
if fid < 0, error('cannot write %s', mdfile); end
c = onCleanup(@() fclose(fid));

w = @(varargin) fprintf(fid, varargin{:});

w('# %s — dependency overview%s', reponame, nl);
w('%s', nl);
if ~isempty(subtitle), w('%s%s%s', subtitle, nl, nl); end

w('> **Generated file — do not edit.** Regenerate with%s', nl);
w('> `%s`%s', regenhint, nl);
w('> (see `clean/LaBGAScore_dep_report.m` in LaBGAScore).%s', nl);
w('> Generated %s by MATLAB %s.%s', ...
    char(datetime('now','Format','yyyy-MM-dd')), version('-release'), nl);
w('%s', nl);

w('This document records which **external** functions each file calls and which%s', nl);
w('repository those live in. Calls that resolve back into this repository, and%s', nl);
w('MathWorks'' own functions, are omitted.%s%s', nl, nl);

if isinf(maxdepth)
    w('Full transitive closure.%s%s', nl, nl);
else
    w('Direct calls only (depth %d). The transitive closure of a second-level%s', maxdepth, nl);
    w('script runs to thousands of files and reduces to "CanlabCore and SPM";%s', nl);
    w('it is what the provenance tooling uses, for a different purpose.%s%s', nl, nl);
end

% --- how to read the caveats ------------------------------------------------

w('## How to read this%s%s', nl, nl);
w('Resolution is static, and MATLAB makes some of it genuinely undecidable.%s', nl);
w('Every edge carries a confidence, and the uncertain ones are reported rather%s', nl);
w('than guessed:%s%s', nl, nl);
w('| Confidence | Meaning |%s', nl);
w('|---|---|%s', nl);
w('| `resolved` | exactly one definition exists for this name |%s', nl);
w('| `ambiguous` | several classes define it — e.g. `threshold` exists in `@atlas`, `@glm_map`, `@image_vector` and `@statistic_image`. Deciding which one runs needs type inference these workspace-chained scripts do not support. All candidates are listed. |%s', nl);
w('| `dotcall` | called as `obj.name(...)`, matched to a class method by name. Could also be a struct field. |%s', nl);
w('| `dynamic` | the file uses `feval`/`eval`/`str2func`, so its real call set cannot be recovered statically |%s', nl);
w('| `unparseable` | the file has a syntax error and could not be walked |%s%s', nl, nl);

% --- summary table ----------------------------------------------------------

w('## Summary%s%s', nl, nl);
w('| Script | Domain | Depends on | Direct calls | Caveats |%s', nl);
w('|---|---|---|---:|---:|%s', nl);

for k = 1:height(R.scripts)
    dep = R.scripts.DependsOn(k);
    if strlength(dep) == 0, dep = "—"; end
    w('| `%s` | %s | %s | %d | %d |%s', R.scripts.Script(k), R.scripts.Domain(k), ...
        dep, R.scripts.NDirect(k), R.scripts.NCaveats(k), nl);
end
w('%s', nl);

% --- repository totals ------------------------------------------------------

if ~isempty(R.edges)
    w('## Dependencies by repository%s%s', nl, nl);
    [cnt, rp] = groupcounts(R.edges.Repo);
    [cnt, ix] = sort(cnt,'descend'); rp = rp(ix);
    w('| Repository | Call edges | Distinct functions |%s', nl);
    w('|---|---:|---:|%s', nl);
    for k = 1:numel(rp)
        nfun = numel(unique(R.edges.Callee(R.edges.Repo == rp(k))));
        w('| %s | %d | %d |%s', rp(k), cnt(k), nfun, nl);
    end
    w('%s', nl);
end

% --- per-script detail ------------------------------------------------------

w('## Per-script detail%s%s', nl, nl);

for k = 1:height(R.scripts)

    rows = R.edges(R.edges.Caller == R.scripts.Script(k),:);

    w('### `%s`%s%s', R.scripts.Script(k), nl, nl);
    w('`%s`%s%s', R.scripts.Path(k), nl, nl);

    if isempty(rows)
        w('No external dependencies.%s%s', nl, nl);
        continue
    end

    urepos = unique(rows.Repo);
    for i = 1:numel(urepos)
        sub = rows(rows.Repo == urepos(i),:);
        w('**%s**%s%s', urepos(i), nl, nl);
        [ucall, ia] = unique(sub.Callee);
        for j = 1:numel(ucall)
            conf = sub.Confidence(ia(j));
            cls = sub.Class(ia(j));
            if strlength(cls) > 0, clsstr = sprintf(' *(@%s)*', cls); else, clsstr = ''; end
            if conf == "resolved"
                w('- `%s`%s%s', ucall(j), clsstr, nl);
            else
                n = sum(sub.Callee == ucall(j));
                w('- `%s`%s — `%s`', ucall(j), clsstr, conf);
                if n > 1, w(', %d candidates', n); end
                w('%s', nl);
            end
        end
        w('%s', nl);
    end

end

end


function local_write_yaml(ymlfile, R, reponame, IDX)
% hand-written rather than round-tripped: the file is generated, flat, and
% consumed by the website's refresh script

nl = newline;
fid = fopen(ymlfile,'w');
if fid < 0, error('cannot write %s', ymlfile); end
c = onCleanup(@() fclose(fid));

w = @(varargin) fprintf(fid, varargin{:});

w('# Dependency data for %s.%s', reponame, nl);
w('# GENERATED by LaBGAScore_dep_report - do not edit by hand.%s', nl);
w('generated: "%s"%s', char(datetime('now','Format','yyyy-MM-dd''T''HH:mm:ss')), nl);
w('matlab: "%s"%s', version('-release'), nl);
w('repo: "%s"%s', reponame, nl);
w('n_indexed_files: %d%s', height(IDX.files), nl);
w('scripts:%s', nl);

for k = 1:height(R.scripts)
    rows = R.edges(R.edges.Caller == R.scripts.Script(k),:);
    w('  - name: "%s"%s', R.scripts.Script(k), nl);
    w('    path: "%s"%s', R.scripts.Path(k), nl);
    w('    domain: "%s"%s', R.scripts.Domain(k), nl);
    w('    n_direct: %d%s', R.scripts.NDirect(k), nl);
    w('    n_caveats: %d%s', R.scripts.NCaveats(k), nl);
    if isempty(rows)
        w('    depends_on: []%s', nl);
        continue
    end
    w('    depends_on:%s', nl);
    urepos = unique(rows.Repo);
    for i = 1:numel(urepos)
        sub = rows(rows.Repo == urepos(i),:);
        q = """";
        w('      - repo: "%s"%s', urepos(i), nl);
        w('        calls: [%s]%s', ...
            strjoin(q + unique(sub.Callee)' + q, ', '), nl);
    end
end

end
