function DEP = LaBGAScore_dep_map(target, varargin)
%
%
% *USAGE*
%
% Maps which dependency files a script or function actually reaches, by
% parsing it with mtree, resolving every called name against the index
% built by LaBGAScore_dep_build_index, and walking the result transitively.
%
% Two consumers:
%
% # LaBGAScore_prov_snapshot / _resolve_retrospective, to answer "do this
%   repo's uncommitted local changes actually affect THIS script?" - which
%   turns a blanket "dirty, state unknown" into a per-script verdict
% # the generated DEPENDENCIES.md documentation
%
% matlab.codetools.requiredFilesAndProducts is NOT used, deliberately: it
% aborts on a syntax error anywhere in the transitive closure, and there
% are 33 such files in the CANlab tree (CanlabCore 10, CanlabPrivate 21,
% canlab_single_trials 1, CANlab_help_examples 1). This function records
% those files as 'unparseable' and carries on.
%
%
% *OPTIONS*
%
% * 'index'         index struct from LaBGAScore_dep_build_index; built on
%                   demand (and cached) if not supplied
% * 'maxdepth'      transitive depth, default Inf (1 = direct calls only)
% * 'followrepos'   only recurse into files whose Repo is in this list,
%                   default: everything except the MATLAB toolbox root
% * 'print'         logical, default true
%
%
% *OUTPUT*
%
% * DEP             struct with fields
%                   .edges  table, one row per call edge: Caller,
%                           CallerFile, Callee, CalleeFile, Repo, Class,
%                           Kind, Depth, Confidence, NCandidates
%                   .files  unique dependency files reached (table)
%                   .repos  unique repos reached, with file counts
%                   .caveats table of edges needing human judgement
%                   .target the file(s) analysed
%
%
% *NOTES ON CONFIDENCE*
%
% * resolved     exactly one candidate definition in the index
% * ambiguous    several classes define the name (threshold, apply_mask,
%                predict ...). All candidates are listed. Deciding which
%                one runs needs type inference that these workspace-chained
%                scripts do not support, so it is REPORTED, NOT GUESSED.
% * dotcall      called as obj.name(...); resolved to a class method by
%                name only. May be a struct field access instead.
% * unparseable  the file has a syntax error and could not be walked
% * dynamic      the file uses feval/eval/str2func/evalin, so its true call
%                set cannot be recovered statically. Flagged per file.
% * toolbox      resolves inside matlabroot, i.e. MathWorks code
%
% Call names are identified as mtree ID nodes minus names that the source
% assigns to (regex over assignment left-hand sides, for-loop variables and
% function declaration arguments). A variable that shadows a function name
% and is never assigned in the same file can still produce a false edge.
%
%
% *SEE ALSO*
%
% LaBGAScore_dep_build_index, LaBGAScore_prov_snapshot
%
% -------------------------------------------------------------------------
%
% modified by: Lukas Van Oudenhove & Claude Opus 5
%
% date:   KU Leuven, August, 2026
%
% -------------------------------------------------------------------------
%
% LaBGAScore_dep_map.m         v1.0
%
% last modified: 2026/08/31
%
%
%% PARSE OPTIONS
% -------------------------------------------------------------------------

p = inputParser;
p.addParameter('index', [], @(x) isempty(x) || isstruct(x));
p.addParameter('maxdepth', Inf, @isnumeric);
p.addParameter('followrepos', {}, @iscell);
p.addParameter('print', true, @islogical);
p.parse(varargin{:});
opt = p.Results;

if isempty(opt.index)
    IDX = LaBGAScore_dep_build_index('print', opt.print);
else
    IDX = opt.index;
end

targets = local_expand_target(target);

if isempty(targets)
    error('No .m files found for target %s', string(target));
end


%% WALK EACH TARGET
% -------------------------------------------------------------------------

caller = strings(0,1); callerfile = strings(0,1); callee = strings(0,1);
calleefile = strings(0,1); repo = strings(0,1); class = strings(0,1);
kind = strings(0,1); depth = zeros(0,1); confidence = strings(0,1);
ncand = zeros(0,1); root = strings(0,1);

for tt = 1:numel(targets)

    thistarget = targets{tt};
    [~, targetname] = fileparts(thistarget);

    seen = containers.Map('KeyType','char','ValueType','logical');
    queue = {thistarget};
    queuedepth = 0;

    while ~isempty(queue)

        thisfile = queue{1};
        thisdepth = queuedepth(1);
        queue(1) = [];
        queuedepth(1) = [];

        if isKey(seen, thisfile), continue, end
        seen(thisfile) = true;

        if thisdepth >= opt.maxdepth, continue, end

        [~, thiscaller] = fileparts(thisfile);

        [names, isdot, parseok, isdynamic, allids] = local_called_names(thisfile);

        if ~parseok
            caller(end+1,1) = string(thiscaller); %#ok<AGROW>
            callerfile(end+1,1) = string(thisfile); %#ok<AGROW>
            callee(end+1,1) = ""; %#ok<AGROW>
            calleefile(end+1,1) = ""; %#ok<AGROW>
            repo(end+1,1) = ""; root(end+1,1) = ""; %#ok<AGROW>
            class(end+1,1) = ""; kind(end+1,1) = ""; %#ok<AGROW>
            depth(end+1,1) = thisdepth; %#ok<AGROW>
            confidence(end+1,1) = "unparseable"; %#ok<AGROW>
            ncand(end+1,1) = 0; %#ok<AGROW>
            continue
        end

        if isdynamic
            caller(end+1,1) = string(thiscaller); %#ok<AGROW>
            callerfile(end+1,1) = string(thisfile); %#ok<AGROW>
            callee(end+1,1) = ""; calleefile(end+1,1) = ""; %#ok<AGROW>
            repo(end+1,1) = ""; root(end+1,1) = ""; %#ok<AGROW>
            class(end+1,1) = ""; kind(end+1,1) = ""; %#ok<AGROW>
            depth(end+1,1) = thisdepth; %#ok<AGROW>
            confidence(end+1,1) = "dynamic"; %#ok<AGROW>
            ncand(end+1,1) = 0; %#ok<AGROW>
        end

        [inplayrepos, unshadowed] = local_inplay_repos(names(~isdot), IDX);

        for n = 1:numel(names)

            thisname = char(names(n));
            key = lower(thisname);
            iscorebuiltin = false;

            if ~isKey(IDX.map, key), continue, end

            rows = IDX.map(key);
            cands = IDX.files(rows,:);

            % the index is keyed lowercase so lookups are cheap, but MATLAB
            % resolves names case-sensitively: Branch must not match
            % @xmltree/branch.m
            cands = cands(cands.Name == string(thisname),:);
            if isempty(cands), continue, end

            % a mex file for another platform cannot be what ran here, and
            % counting all four builds of liblinear's predict as separate
            % candidates would inflate the ambiguity for no reason
            ismex = cands.Kind == "mex";
            if any(ismex)
                wrongplatform = ismex & ~endsWith(cands.File, ['.' mexext]);
                cands = cands(~wrongplatform,:);
                if isempty(cands), continue, end
            end

            % a dot-call can only sensibly be a class method
            if isdot(n)
                cands = cands(cands.Kind == "method",:);
                if isempty(cands), continue, end

                % obj.name is indistinguishable from a property or struct
                % field read, and object PROPERTIES are the common case in
                % this code: atlas_obj.labels, V.dim, V.fname. Matching those
                % by name alone invented dependencies on BrainSpace (labels),
                % gift (dim, fname) and MKDA (test), none of which these
                % scripts use. A method call cannot introduce a repository
                % that nothing else in the file reaches - you must construct
                % or load the object first, and that is a plain call - so a
                % dot-call may only reinforce a repo already in play.
                cands = cands(ismember(cands.Repo, inplayrepos),:);
                if isempty(cands), continue, end
            elseif any(startsWith(cands.File, fullfile(matlabroot,'toolbox','matlab')))
                % A core-language function of this name exists (isempty, find,
                % table, split, error - everything under toolbox/matlab). A
                % plain call is that, full stop: CanlabCore defining a
                % @region/table method does not make every table() call a
                % CanlabCore dependency.
                cands = cands(startsWith(cands.File, matlabroot),:);
                if isempty(cands), continue, end
                iscorebuiltin = true;

            elseif numel(unique(cands.Repo)) > 1
                % Genuinely ambiguous ACROSS repositories. Several copies of a
                % name inside one repository (spm12 ships two spm_vol.m) are
                % not ambiguous in the sense that matters here.
                %
                % MATLAB's own path order decides which one actually runs, so
                % ask it. That is exactly the question, and unlike using
                % which() as the primary resolver it cannot miss a class
                % method: the candidate list is already known, which() only
                % picks among it.
                w = local_which_cached(thisname);
                pick = strings(0,1);
                if ~isempty(w)
                    match = cands.File == string(w);
                    if any(match), pick = cands.Repo(find(match,1)); end
                end
                if ~isempty(pick)
                    keep = cands.Repo == pick | ...
                           (cands.Kind == "method" & ismember(cands.Class, allids));
                else
                    keep = startsWith(cands.File, matlabroot) | ...
                           ismember(cands.Repo, inplayrepos) | ...
                           (cands.Kind == "method" & ismember(cands.Class, allids));
                end
                cands = cands(keep,:);
                if isempty(cands), continue, end
                % An ambiguous name is weak evidence: several repos define it
                % and we cannot tell which runs. Let it reinforce a repo the
                % file already reaches, but never introduce a new one -
                % otherwise spm_defaults pulls in cocoanCORE, get_var pulls in
                % ooFmriDataObjML, and so on for toolboxes nothing here uses.
                keep = startsWith(cands.File, matlabroot) | ...
                       ismember(cands.Repo, inplayrepos) | ...
                       (cands.Kind == "method" & ismember(cands.Class, allids));
                cands = cands(keep,:);
                if isempty(cands), continue, end

            elseif any(startsWith(cands.File, matlabroot)) && ~ismember(key, unshadowed)
                % A MathWorks function of this name exists, so a plain call is
                % most likely that. Keep a third-party candidate only when it
                % is a method of a class the file actually references.
                % Without this, every table(), find() or error() picks up
                % CanlabCore's @region/table.m and friends and invents a
                % dependency on half the CANlab tree. The test is deliberately
                % narrow: names with no builtin of their own - merge_atlases,
                % select_atlas_subset - are left alone entirely.
                keep = startsWith(cands.File, matlabroot) | ...
                       (cands.Kind == "method" & ismember(cands.Class, allids)) | ...
                       ismember(cands.Repo, inplayrepos);
                cands = cands(keep,:);
                if isempty(cands), continue, end
            end

            % never treat the file as calling itself
            cands = cands(cands.File ~= string(thisfile),:);
            if isempty(cands), continue, end

            istoolbox = startsWith(cands.File, matlabroot);

            if iscorebuiltin
                thisconf = "toolbox";
            elseif isdot(n)
                thisconf = "dotcall";
            elseif numel(unique(cands.Repo)) > 1
                thisconf = "ambiguous";
            elseif height(cands) > 1
                thisconf = "ambiguous_within_repo";
            elseif all(istoolbox)
                thisconf = "toolbox";
            else
                thisconf = "resolved";
            end

            for c = 1:height(cands)
                caller(end+1,1) = string(thiscaller); %#ok<AGROW>
                callerfile(end+1,1) = string(thisfile); %#ok<AGROW>
                callee(end+1,1) = string(thisname); %#ok<AGROW>
                calleefile(end+1,1) = cands.File(c); %#ok<AGROW>
                repo(end+1,1) = cands.Repo(c); %#ok<AGROW>
                root(end+1,1) = cands.Root(c); %#ok<AGROW>
                class(end+1,1) = cands.Class(c); %#ok<AGROW>
                kind(end+1,1) = cands.Kind(c); %#ok<AGROW>
                depth(end+1,1) = thisdepth + 1; %#ok<AGROW>
                if istoolbox(c)
                    confidence(end+1,1) = "toolbox"; %#ok<AGROW>
                else
                    confidence(end+1,1) = thisconf; %#ok<AGROW>
                end
                ncand(end+1,1) = height(cands); %#ok<AGROW>
            end

            % queue for the transitive walk
            for c = 1:height(cands)
                if istoolbox(c), continue, end
                if ~endsWith(cands.File(c),'.m'), continue, end
                if ~isempty(opt.followrepos) && ~ismember(char(cands.Repo(c)), opt.followrepos)
                    continue
                end
                if ~isKey(seen, char(cands.File(c)))
                    queue{end+1} = char(cands.File(c)); %#ok<AGROW>
                    queuedepth(end+1) = thisdepth + 1; %#ok<AGROW>
                end
            end

        end

    end

    if opt.print
        fprintf('  %-56s %4d edges\n', targetname, sum(callerfile ~= "" & strcmp(cellstr(caller), caller)));
    end

end

DEP.edges = table(caller, callerfile, callee, calleefile, repo, root, class, ...
    kind, depth, confidence, ncand, ...
    'VariableNames', {'Caller','CallerFile','Callee','CalleeFile','Repo', ...
                      'Root','Class','Kind','Depth','Confidence','NCandidates'});

% a name can be picked up both as a plain identifier and as a dot-call,
% which would otherwise emit the same edge twice
[~, keep] = unique(DEP.edges(:,{'CallerFile','Callee','CalleeFile'}), 'rows', 'stable');
DEP.edges = DEP.edges(keep,:);

DEP.target = targets;


%% SUMMARISE
% -------------------------------------------------------------------------

real = DEP.edges(DEP.edges.CalleeFile ~= "",:);

[uf, ia] = unique(real.CalleeFile);
DEP.files = table(uf, real.Repo(ia), real.Class(ia), real.Kind(ia), ...
    'VariableNames', {'File','Repo','Class','Kind'});

if isempty(uf)
    DEP.repos = table(strings(0,1), zeros(0,1), 'VariableNames', {'Repo','NFiles'});
else
    [ur, ~, ix] = unique(DEP.files.Repo);
    DEP.repos = table(ur, accumarray(ix,1), 'VariableNames', {'Repo','NFiles'});
    DEP.repos = sortrows(DEP.repos,'NFiles','descend');
end

DEP.caveats = DEP.edges(ismember(DEP.edges.Confidence, ...
    ["ambiguous","dotcall","unparseable","dynamic"]),:);


%% PRINT REPORT
% -------------------------------------------------------------------------

if opt.print
    fprintf('\n%d file(s) analysed, %d edges, %d distinct dependency files\n', ...
        numel(targets), height(DEP.edges), height(DEP.files));
    fprintf('\nRepositories reached:\n');
    for k = 1:height(DEP.repos)
        fprintf('  %-32s %5d files\n', DEP.repos.Repo(k), DEP.repos.NFiles(k));
    end
    if ~isempty(DEP.caveats)
        fprintf('\n%d edge(s) need human judgement:\n', height(DEP.caveats));
        for c = ["ambiguous","dotcall","unparseable","dynamic"]
            n = sum(DEP.caveats.Confidence == c);
            if n > 0, fprintf('  %-14s %d\n', c, n); end
        end
        fprintf('See DEP.caveats. These are reported, not guessed.\n');
    end
    fprintf('\n');
end

end


%% LOCAL FUNCTIONS
% -------------------------------------------------------------------------

function files = local_expand_target(target)
% accept a file, a folder, or a cellstr/string array of either

if ischar(target) || (isstring(target) && isscalar(target))
    target = {char(target)};
else
    target = cellstr(target);
end

files = {};

for k = 1:numel(target)
    t = target{k};
    if isfolder(t)
        d = dir(fullfile(t,'**','*.m'));
        d = d(~[d.isdir]);
        for i = 1:numel(d)
            files{end+1} = fullfile(d(i).folder, d(i).name); %#ok<AGROW>
        end
    elseif isfile(t)
        files{end+1} = t; %#ok<AGROW>
    end
end

end


function [names, isdot, parseok, isdynamic, allids] = local_called_names(fullpath)
% mtree ID and FIELD nodes, minus anything the file assigns to

names = strings(0,1); isdot = false(0,1); parseok = false; isdynamic = false;
allids = strings(0,1);

try
    t = mtree(fullpath, '-file');
catch
    return
end

if t.anykind('ERR')
    return
end

parseok = true;

% stringvals returns a cellstr whose orientation varies, so force columns
ids = string(mtfind(t,'Kind','ID').stringvals());
fields = string(mtfind(t,'Kind','FIELD').stringvals());

ids = ids(:);
fields = fields(:);

if isempty(ids), ids = strings(0,1); end
if isempty(fields), fields = strings(0,1); end

isdynamic = any(ismember(ids, ["feval","eval","evalin","str2func","assignin","evalc"]));

assigned = local_assigned_names(fullpath);

ids = unique(ids);
ids = ids(~ismember(ids, assigned));
ids = ids(:);

fields = unique(fields);
fields = fields(:);

% obj.name(...) is only a method call if obj is a MATLAB object. Two common
% cases where it is not, and where matching by name alone invents
% dependencies: struct fields, and Java objects.
[structfields, javafields] = local_nonmethod_fields(fullpath);
fields = fields(~ismember(fields, [structfields; javafields; assigned]));

% In a classdef, the names declared in a properties block are identifiers as
% far as mtree is concerned, but they are declarations, not calls. Without
% this, CanlabCore's @atlas property "labels" resolved to BrainSpace's
% labels.m and invented a dependency on a toolbox nothing here uses.
declared = local_classdef_declarations(fullpath);
ids = ids(~ismember(ids, declared));

names = [ids; fields];
isdot = [false(numel(ids),1); true(numel(fields),1)];

% every identifier and field in the file, used to decide whether a class is
% in play at all (see the method filter in the caller)
allids = unique([ids; fields]);

end


function w = local_which_cached(name)
% which() is a path lookup per call; the same names recur constantly across
% a transitive walk, so remember the answers.

persistent cache
if isempty(cache), cache = containers.Map('KeyType','char','ValueType','any'); end

if isKey(cache, name)
    w = cache(name);
    return
end

try
    w = which(name);
catch
    w = '';
end

cache(name) = w;

end


function [repos, unshadowed] = local_inplay_repos(names, IDX)
% Repositories this file demonstrably uses, judged only on PLAIN calls whose
% name no MathWorks function shadows - load_atlas, merge_atlases, fmri_data
% and so on. Dot-calls are excluded deliberately: they are gated on this
% result, so letting them contribute would make the test circular.
% Those are unambiguous evidence. Once a repository is established as in play,
% a shadowed name like predict or threshold can reasonably be attributed to it
% as well; for a file that never touches the repository at all, it cannot.

repos = strings(0,1);
unshadowed = strings(0,1);

for n = 1:numel(names)
    key = lower(char(names(n)));
    if ~isKey(IDX.map, key), continue, end
    c = IDX.files(IDX.map(key),:);
    c = c(c.Name == names(n),:);
    if isempty(c), continue, end
    if any(startsWith(c.File, matlabroot)), continue, end
    % Evidence only if the name belongs to exactly ONE repository. A name
    % several repos define - spm_defaults lives in spm12, cocoanCORE AND
    % CanlabPrivate - says nothing about which of them this file uses, and
    % treating it as evidence put all three in play, which is how cocoanCORE
    % entered records for scripts that never touch it.
    if numel(unique(c.Repo)) > 1, continue, end
    unshadowed(end+1,1) = string(key); %#ok<AGROW>
    repos = [repos; c.Repo]; %#ok<AGROW>
end

repos = unique(repos);
unshadowed = unique(unshadowed);

end


function declared = local_classdef_declarations(fullpath)
% Names declared in a classdef's properties (and events) blocks. mtree
% reports them as ID nodes, but they are declarations rather than calls.

declared = strings(0,1);

src = string(fileread(fullpath));
lines = splitlines(src);

if ~any(~cellfun(@isempty, regexp(cellstr(strtrim(lines)), '^classdef\>', 'once')))
    return                                  % not a classdef, nothing to do
end

depth = 0;
inblock = false;

for k = 1:numel(lines)

    raw = strtrim(lines(k));
    if startsWith(raw, '%') || strlength(raw) == 0, continue, end

    % strip trailing comment
    line = regexprep(raw, '%.*$', '');
    line = strtrim(line);
    if strlength(line) == 0, continue, end

    if ~inblock
        if ~isempty(regexp(line, '^(properties|events)\>', 'once'))
            inblock = true; depth = 1;
        end
        continue
    end

    if ~isempty(regexp(line, '^(if|for|while|switch|try|function|parfor)\>', 'once'))
        depth = depth + 1;
    elseif ~isempty(regexp(line, '^end\>', 'once'))
        depth = depth - 1;
        if depth <= 0, inblock = false; end
        continue
    end

    nm = regexp(line, '^([A-Za-z]\w*)', 'tokens', 'once');
    if ~isempty(nm)
        declared(end+1,1) = string(nm{1}); %#ok<AGROW>
    end

end

declared = unique(declared);

end


function [structfields, javafields] = local_nonmethod_fields(fullpath)
% field names that cannot be class-method calls

src = string(fileread(fullpath));
lines = splitlines(src);
lines = lines(~startsWith(strtrim(lines),'%'));
txt = strjoin(lines, newline);

% anything ever assigned as x.field = ... is a struct field
structfields = strings(0,1);
tok = regexp(txt, '[A-Za-z]\w*(?:\([^)]*\))?\.([A-Za-z]\w*)\s*=(?![=])', 'tokens');
for k = 1:numel(tok)
    structfields(end+1,1) = string(tok{k}{1}); %#ok<AGROW>
end
structfields = unique(structfields);

% variables assigned from a Java constructor: their dot-calls are Java
% methods (md.update, jf.length, ...), not MATLAB class methods
javafields = strings(0,1);
tok = regexp(txt, '([A-Za-z]\w*)\s*=\s*(?:java(?:Object)?[\.\(])', 'tokens');
javavars = strings(0,1);
for k = 1:numel(tok)
    javavars(end+1,1) = string(tok{k}{1}); %#ok<AGROW>
end
javavars = unique(javavars);
for k = 1:numel(javavars)
    t2 = regexp(txt, char(javavars(k)) + "\.([A-Za-z]\w*)", 'tokens');
    for i = 1:numel(t2)
        javafields(end+1,1) = string(t2{i}{1}); %#ok<AGROW>
    end
end
javafields = unique(javafields);

end


function assigned = local_assigned_names(fullpath)
% names the file assigns to: LHS of =, for-loop variables, function args.
% Regex over the source rather than an mtree subtree walk, which mtree does
% not expose conveniently. Over-inclusive by design: dropping a genuine call
% because a variable shares its name is preferable to inventing a
% dependency on a repo the script never touches.

assigned = strings(0,1);

src = fileread(fullpath);
lines = splitlines(string(src));

% strip full-line comments; keeps it simple and is enough for LHS detection
lines = lines(~startsWith(strtrim(lines),'%'));

% simple assignment:  name = ...   /  name(idx) = ...  /  name.field = ...
tok = regexp(lines, '^\s*([A-Za-z]\w*)\s*[\(\{\.]?[^=]*=(?![=])', 'tokens');
for k = 1:numel(tok)
    if ~isempty(tok{k}), assigned(end+1,1) = string(tok{k}{1}); end %#ok<AGROW>
end

% multiple assignment:  [a, b, ~] = ...
tok = regexp(lines, '^\s*\[([^\]]*)\]\s*=(?![=])', 'tokens');
for k = 1:numel(tok)
    if isempty(tok{k}), continue, end
    parts = strsplit(tok{k}{1}, {',',' '});
    for i = 1:numel(parts)
        nm = regexp(parts{i}, '^([A-Za-z]\w*)', 'tokens', 'once');
        if ~isempty(nm), assigned(end+1,1) = string(nm{1}); end %#ok<AGROW>
    end
end

% for-loop variables
tok = regexp(lines, '^\s*(?:for|parfor)\s*\(?\s*([A-Za-z]\w*)\s*=', 'tokens');
for k = 1:numel(tok)
    if ~isempty(tok{k}), assigned(end+1,1) = string(tok{k}{1}); end %#ok<AGROW>
end

% function declaration: outputs and inputs
tok = regexp(lines, '^\s*function\s+(.*)$', 'tokens');
for k = 1:numel(tok)
    if isempty(tok{k}), continue, end
    decl = tok{k}{1};
    nm = regexp(decl, '[A-Za-z]\w*', 'match');
    % first name after any "=" is the function itself; keep all as assigned
    for i = 1:numel(nm)
        assigned(end+1,1) = string(nm{i}); %#ok<AGROW>
    end
end

assigned = unique(assigned);

end
