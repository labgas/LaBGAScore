function [T, G] = LaBGAScore_firstlevel_task_motion_diagnostics(BIDSdir, derivdir, varargin)
% Quantify task-correlated head motion and its leakage into first-level betas
%
% :Usage:
% ::
%
%     [T, G] = LaBGAScore_firstlevel_task_motion_diagnostics(BIDSdir, derivdir, 'conditions', C, ...)
%
% Answers two questions about a first-level design, without refitting anything
% and without touching the imaging data:
%
%   1. Is head motion correlated with the task, and is that correlation
%      consistent in sign across subjects? Only a consistent one survives
%      group averaging and can bias a group map; an inconsistent one adds
%      between-subject variance and costs power instead.
%
%   2. How far can a motion artifact displace each contrast estimate, in
%      units of that contrast's own standard error, given the nuisance
%      regressors actually in the model?
%
% Both are computed from the fMRIPrep confound files and the BIDS events
% files only, so it runs in seconds and needs no SPM.mat and no betas.
%
% :Inputs:
%
%   **BIDSdir:**
%        BIDS directory. Events files are resolved per run by the BIDS
%        inheritance principle, so both layouts work without an option: one
%        events file per subject and run under sub-*/func/, or a single
%        task-level file such as BIDSdir/task-<label>_events.tsv shared by
%        every run, which is how a long block design with fixed timing is
%        usually stored. The most specific applicable file wins, and the route
%        taken is printed once per run of the function.
%
%   **derivdir:**
%        fMRIPrep directory holding sub-*/func/*_desc-confounds_timeseries.tsv
%
% :Optional Inputs:
%
%   **'conditions':**
%        cell array of trial_type values to model. Default: every unique
%        trial_type found across the events files
%
%   **'dsgn':**
%        the CANlab DSGN struct, as built by
%        <study>_firstlevel_m<M>_s1_options_dsgn_struct. Conditions are taken
%        from DSGN.conditions{1} and contrasts from DSGN.contrasts,
%        DSGN.contrastweights and DSGN.contrastnames, so neither needs to be
%        written out by hand. Note that DSGN.contrastweights cannot simply be
%        flattened with cell2mat: its vectors differ in length and its entries
%        are weights on regexp GROUPS rather than on conditions, so they are
%        mapped through DSGN.contrasts. The derived matrix is printed for
%        checking, and contrasts that cannot be mapped - those on parametric
%        modulators, which this function does not model - are dropped with a
%        warning. 'conditions' and 'contrasts' override it if both are given
%
%   **'contrasts':**
%        [nContrast x nCondition] weight matrix. Default: from 'dsgn' if given,
%        otherwise one contrast per condition against implicit baseline
%
%   **'contrastnames':**
%        cell array of names for the rows of 'contrasts'
%
%   **'tr':**
%        repetition time in seconds. REQUIRED unless it can be read from a
%        sidecar .json next to the events file
%
%   **'hpf':**
%        high-pass filter length in seconds, default 180 (the LaBGAS default;
%        must match DSGN.hpf, since SPM filters the design matrix too)
%
%   **'quadratic':**
%        true if the model included quadratic motion terms (24 columns),
%        false for the 12-column version. Default true
%
%   **'task':**
%        task label to select files with, e.g. 'sweettaste'. Default: all
%
%   **'subjects':**
%        cell array of subject directory names. Default: all in derivdir
%
% :Outputs:
%
%   **T:**
%        table, one row per subject x run x condition, with the correlation
%        between the convolved regressor and framewise displacement and each
%        of the six realignment parameters, and the exposure metric
%
%   **G:**
%        struct of group-level summaries: one-sample t across subjects on the
%        Fisher-z correlations (the consistency test), and the per-contrast
%        leakage bound
%
% :Notes:
%
% The leakage bound is the fraction of the contrast estimator's weight vector
% that a motion artifact can act through, given the nuisance columns that were
% actually fitted. Multiply it by the artifact-to-noise ratio to read it as a
% shift in standard errors. It is a worst case over artifact directions: a
% single artifact direction shared between two conditions cancels in their
% difference, and the bound does not credit that cancellation.
%
% :See also: LaBGAScore_firstlevel_s2_fit_model, scn_spm_design_check
%
% -------------------------------------------------------------------------
%
% authors: Lukas Van Oudenhove
%
% date: KU Leuven, September, 2026
%
% -------------------------------------------------------------------------
%
% LaBGAScore_firstlevel_task_motion_diagnostics.m         v1.3
%
% last modified: 2026/09/03


%% PARSE OPTIONS
% -------------------------------------------------------------------------

p = inputParser;
p.addParameter('conditions', {}, @iscell);
p.addParameter('contrasts', [], @isnumeric);
p.addParameter('contrastnames', {}, @iscell);
p.addParameter('tr', [], @isnumeric);
p.addParameter('hpf', 180, @isnumeric);
p.addParameter('quadratic', true, @islogical);
p.addParameter('task', '', @(x) ischar(x) || isstring(x));
p.addParameter('subjects', {}, @iscell);
p.addParameter('dsgn', [], @(x) isempty(x) || isstruct(x));
p.parse(varargin{:});
opt = p.Results;

if isempty(opt.tr), error('tr is required: pass ''tr'', DSGN.tr'); end

subdirs = dir(fullfile(derivdir, 'sub-*'));
subdirs = subdirs([subdirs.isdir]);
subs = {subdirs.name}';
if ~isempty(opt.subjects), subs = intersect(subs, opt.subjects, 'stable'); end
if isempty(subs), error('no sub-* directories found in %s', derivdir); end

taskpat = char(opt.task); if ~isempty(taskpat), taskpat = ['*' taskpat '*']; else, taskpat = '*'; end


%% COLLECT CONDITIONS IF NOT GIVEN
% -------------------------------------------------------------------------

conds = opt.conditions;
if isempty(conds) && ~isempty(opt.dsgn) && isfield(opt.dsgn, 'conditions') ...
        && ~isempty(opt.dsgn.conditions)
    % DSGN.conditions holds one cell per session, and in a multisession design
    % those differ: the drug session has one set of conditions and the placebo
    % session another. Taking only the first would silently drop half the
    % design and every contrast that refers to it, so take the union, keeping
    % the order in which conditions first appear.
    conds = {};
    for i = 1:numel(opt.dsgn.conditions)
        ci = opt.dsgn.conditions{i};
        if ~iscell(ci), ci = {ci}; end
        conds = [conds, setdiff(ci, conds, 'stable')]; %#ok<AGROW>
    end
end
if isempty(conds)
    % look wherever events files are allowed to live, not only under the
    % subject directories, so a single task-level file is picked up too
    e = [ dir(fullfile(BIDSdir, [taskpat '_events.tsv'])); ...
          dir(fullfile(BIDSdir, 'sub-*', [taskpat '_events.tsv'])); ...
          dir(fullfile(BIDSdir, 'sub-*', 'func', [taskpat '_events.tsv'])); ...
          dir(fullfile(BIDSdir, 'sub-*', 'ses-*', 'func', [taskpat '_events.tsv'])) ];
    e = e(~[e.isdir]' & ~contains({e.name}, 'noninterest')');
    seen = string([]);
    for f = 1:numel(e)
        ev = readtable(fullfile(e(f).folder, e(f).name), 'FileType','text', ...
            'Delimiter','\t', 'VariableNamingRule','preserve');
        if ~ismember('trial_type', ev.Properties.VariableNames), continue, end
        seen = union(seen, unique(string(ev.trial_type)));
    end
    if isempty(seen)
        error('no events files with a trial_type column found under %s', BIDSdir);
    end
    conds = cellstr(seen);
end
nC = numel(conds);

CW = opt.contrasts; cn = opt.contrastnames;
if isempty(CW) && ~isempty(opt.dsgn)
    [CW, cn] = local_contrasts_from_dsgn(opt.dsgn, conds);
end
if isempty(CW), CW = eye(nC); end
if isempty(cn), cn = arrayfun(@(k) sprintf('contrast %d', k), 1:size(CW,1), 'UniformOutput', false); end


%% LOOP OVER SUBJECTS AND RUNS
% -------------------------------------------------------------------------

rows = {}; Lk = []; Lsub = {}; Lidx = []; evreported = {};
mot = {'trans_x','trans_y','trans_z','rot_x','rot_y','rot_z'};

for s = 1:numel(subs)

    % multisession studies put the runs under sub-*/ses-*/func rather than
    % sub-*/func, so look in both rather than silently finding nothing
    cfs = [ dir(fullfile(derivdir, subs{s}, 'func', [taskpat '_desc-confounds_timeseries.tsv'])); ...
            dir(fullfile(derivdir, subs{s}, 'ses-*', 'func', [taskpat '_desc-confounds_timeseries.tsv'])) ];

    for f = 1:numel(cfs)

        stem = erase(cfs(f).name, '_desc-confounds_timeseries.tsv');
        [evf, how] = local_find_events(BIDSdir, subs{s}, stem);
        if isempty(evf)
            warning('no events file resolved for %s, skipping', stem); continue
        end
        if ~strcmp(how, 'exact') && ~any(strcmp(how, evreported))
            fprintf('events files resolved by %s (e.g. %s -> %s)\n', how, stem, evf);
            evreported{end+1} = how; %#ok<AGROW>
        end

        cf = readtable(fullfile(cfs(f).folder, cfs(f).name), 'FileType','text', ...
            'Delimiter','\t', 'TreatAsEmpty','n/a', 'VariableNamingRule','preserve');
        ev = readtable(evf, 'FileType','text', 'Delimiter','\t', 'VariableNamingRule','preserve');

        n  = height(cf);
        vn = cf.Properties.VariableNames;
        if ~all(ismember(mot, vn))
            warning('%s lacks the six realignment parameters, skipping', stem); continue
        end

        R = cf{:, mot}; R(~isfinite(R)) = 0;
        fd = zeros(n,1);
        if ismember('framewise_displacement', vn), fd = cf.framewise_displacement; fd(~isfinite(fd)) = 0; end
        csf = zeros(n,0);
        if ismember('csf', vn), csf = cf.csf; csf(~isfinite(csf)) = 0; end
        sp = table2array(cf(:, contains(vn, 'motion_outlier')));

        d1 = local_grad(R);
        d2 = local_grad(d1);
        zs = @(X) (X - mean(X)) ./ max(std(X), eps);

        A = zs([R d1]);            % the intended Friston-style expansion
        B = [d1 d2];               % what LaBGAScore fitted before the 2026/09/02 fix
        if opt.quadratic, A = [A A.^2]; B = [B B.^2]; end %#ok<AGROW>

        X = local_convolve(ev, conds, n, opt.tr);
        K = local_dct(n, opt.tr, opt.hpf);

        QA = orth(A - mean(A));
        Df = [X K sp csf B ones(n,1)];      % the full design, task columns included
        Dpf = pinv(Df);

        % Exposure is asked of the NUISANCE part of the design only. Including
        % the task columns here would project the task regressor onto a space
        % that already contains it, and every exposure would come back zero.
        Fn = [K sp csf B ones(n,1)];
        QN = orth(Fn);
        M  = A - QN*(QN'*A); M = M(:, sum(abs(M)) > 1e-10); QM = orth(M);

        % In a multisession design the conditions of one session are absent from
        % another, leaving all-zero columns here. Those are not estimable in
        % this run and neither is any contrast that loads on them.
        present = std(X) > eps;

        Kr = @(v) v - K*(K\v);
        for c = 1:nC
            if ~present(c), continue, end
            x = X(:,c) - mean(X(:,c));
            xr = x - QN*(QN'*x);            % task variance surviving the fitted nuisance model
            rows(end+1,:) = { subs{s}, stem, conds{c}, ...
                norm(QM'*xr)^2 / max(norm(xr)^2, eps), ...
                corr(Kr(x), Kr(fd)), num2cell(corr(Kr(x), Kr(R)))  }; %#ok<AGROW>
        end

        % per-contrast leakage bound
        for k = 1:size(CW,1)
            if any(CW(k, ~present) ~= 0), continue, end   % spans a session this run lacks
            w = zeros(1, size(Df,2)); w(1:nC) = CW(k,:);
            r = w * Dpf;
            if norm(r) < eps, continue, end
            Lk(end+1,1) = norm(r*QA) / norm(r); %#ok<AGROW>
            Lsub{end+1,1} = subs{s}; Lidx(end+1,1) = k; %#ok<AGROW>
        end
    end
end

if isempty(rows)
    error(['no runs could be analysed. Looked for confound files matching\n' ...
           '  %s\nand\n  %s\nCheck that derivdir points at the fmriprep directory ' ...
           'holding sub-*, and that the ''task'' label matches the filenames.'], ...
           fullfile(derivdir, 'sub-*', 'func', [taskpat '_desc-confounds_timeseries.tsv']), ...
           fullfile(derivdir, 'sub-*', 'ses-*', 'func', [taskpat '_desc-confounds_timeseries.tsv']));
end

rmot = cell2mat(cellfun(@cell2mat, rows(:,6), 'UniformOutput', false));
T = table(rows(:,1), rows(:,2), rows(:,3), cell2mat(rows(:,4)), cell2mat(rows(:,5)), rmot, ...
    'VariableNames', {'subject','run','condition','exposure','r_fd','r_motion'});


%% GROUP SUMMARIES
% -------------------------------------------------------------------------

usub = unique(T.subject, 'stable'); nS = numel(usub);
G = struct();
G.conditions = conds(:)';
G.consistency = table();

fz = @(r) atanh(min(max(r, -0.999), 0.999));
tstat = nan(nC,1); pval = nan(nC,1); rbar = nan(nC,1); expo = nan(nC,1);
for c = 1:nC
    z = nan(nS,1); e = nan(nS,1);
    for s = 1:nS
        m = strcmp(T.subject, usub{s}) & strcmp(T.condition, conds{c});
        if ~any(m), continue; end
        z(s) = mean(fz(T.r_fd(m))); e(s) = mean(T.exposure(m));
    end
    z = z(~isnan(z)); e = e(~isnan(e));
    if numel(z) > 2, [~, pval(c), ~, st] = ttest(z); tstat(c) = st.tstat; end
    rbar(c) = tanh(mean(z)); expo(c) = mean(e);
end
G.consistency = table(conds(:), rbar, tstat, pval, expo, ...
    'VariableNames', {'condition','mean_r_fd','t_across_subjects','p','mean_exposure'});

lk = nan(size(CW,1), 1); lkmax = nan(size(CW,1), 1);
for k = 1:size(CW,1)
    per = nan(nS,1);
    for s = 1:nS, per(s) = mean(Lk(strcmp(Lsub, usub{s}) & Lidx == k)); end
    lk(k) = mean(per, 'omitnan'); lkmax(k) = max(per, [], 'omitnan');
end
G.leakage = table(cn(:), lk, lkmax, 'VariableNames', {'contrast','mean_leakage_SE','max_subject'});

fprintf('\n%d subjects, %d runs\n\n', nS, numel(unique(T.run)));
disp(G.consistency);
fprintf(['\nmean_r_fd is the correlation between the convolved regressor and framewise\n' ...
         'displacement; a large |t_across_subjects| means it points the same way in most\n' ...
         'subjects and so survives group averaging.\n\n']);
disp(G.leakage);
fprintf(['\nmean_leakage_SE: shift in the contrast estimate, in units of its own standard\n' ...
         'error, per unit ratio of motion artifact to residual noise. Zero if the nuisance\n' ...
         'model spans the artifact.\n']);

nanrows = find(isnan(G.leakage.mean_leakage_SE));
if ~isempty(nanrows)
    fprintf(['\n%d contrast(s) came back NaN: %s.\n' ...
             'This diagnostic works one run at a time, so a contrast comparing conditions\n' ...
             'that live in different sessions is not estimable within any single run. That\n' ...
             'is a limit of the method, not a property of the design; such contrasts have to\n' ...
             'be judged from the per-condition rows above.\n'], ...
             numel(nanrows), strjoin(G.leakage.contrast(nanrows)', ', '));
end

end % main


%% LOCAL FUNCTIONS
% -------------------------------------------------------------------------

function d = local_grad(X)
    d = cell2mat(arrayfun(@(c) gradient(X(:,c)), 1:size(X,2), 'UniformOutput', false));
end

function idx = local_match_conditions(pat, conds)
% Which conditions does one of DSGN's contrast regexps select?
%
% The patterns are written to match beta regressor names, not bare condition
% names: '.*stress{1}\s[^x]' wants "stress", one whitespace, then a character
% that is not x, which is how the unmodulated regressor is told apart from the
% "stress x<pmod>" interaction. A bare condition name has no trailing
% whitespace and so never matches. Each condition is therefore tested in
% several renderings of how SPM might name its regressor.
    idx = [];
    for k = 1:numel(conds)
        c = conds{k};
        renders = { c, [c ' '], [c ' *bf(1)'], ...
                    ['Sn(1) ' c '*bf(1)'], ['Sn(1) ' c ' *bf(1)'] };
        for q = 1:numel(renders)
            if ~isempty(regexp(renders{q}, pat, 'once'))
                idx(end+1) = k; %#ok<AGROW>
                break
            end
        end
    end
end

function [CW, cn] = local_contrasts_from_dsgn(DSGN, conds)
% Build a [nContrast x nCondition] weight matrix from a CANlab DSGN struct.
%
% DSGN.contrastweights is a cell array of numeric vectors of DIFFERENT lengths,
% so cell2mat cannot flatten it, and its entries are not per condition in any
% case. Each weight applies to one regexp GROUP in DSGN.contrasts{c}, which
% selects a set of beta regressors:
%
%   DSGN.contrasts{c}       = {{'.*stress{1}\s[^x]'} {'.*control{1}\s[^x]'}}
%   DSGN.contrastweights{c} = [1 -1]
%
% so weight 1 goes to whatever the first pattern matches and -1 to the second.
% The weight is assigned here to every condition that group's pattern selects.
%
% Contrasts on regressors this function does not model - parametric modulators
% above all, since it builds regressors from trial_type alone - match no
% condition and are dropped with a warning rather than silently mis-mapped.
    if ~isfield(DSGN, 'contrasts') || ~isfield(DSGN, 'contrastweights')
        error('dsgn needs both .contrasts and .contrastweights');
    end
    nCon = numel(DSGN.contrasts); nC = numel(conds);
    CW = zeros(nCon, nC); cn = cell(nCon, 1); keep = true(nCon, 1); why = cell(nCon,1);

    for c = 1:nCon
        if isfield(DSGN, 'contrastnames') && numel(DSGN.contrastnames) >= c
            cn{c} = DSGN.contrastnames{c};
        else
            cn{c} = sprintf('contrast %d', c);
        end

        grp = DSGN.contrasts{c};
        w   = DSGN.contrastweights{c};
        if iscell(w), w = cell2mat(w); end
        if ~isnumeric(w)
            keep(c) = false; why{c} = 'contrastweights is not numeric'; continue
        end
        if numel(w) ~= numel(grp)
            keep(c) = false;
            why{c} = sprintf('%d weights for %d regexp groups', numel(w), numel(grp));
            continue
        end

        for j = 1:numel(grp)
            pats = grp{j};
            if ~iscell(pats), pats = {pats}; end
            idx = [];
            for q = 1:numel(pats)
                idx = union(idx, local_match_conditions(pats{q}, conds));
            end
            if isempty(idx)
                keep(c) = false;
                why{c} = sprintf('pattern "%s" matches none of the conditions', pats{1});
                break
            end
            CW(c, idx) = CW(c, idx) + w(j);
        end
    end

    if any(~keep)
        msg = '';
        for c = find(~keep)'
            msg = [msg sprintf('\n    %-32s %s', cn{c}, why{c})]; %#ok<AGROW>
        end
        warning(['%d of %d contrasts in dsgn could not be mapped onto the conditions ' ...
                 'and are excluded:%s\nPass ''contrasts'' explicitly to override.'], ...
                 sum(~keep), nCon, msg);
    end
    CW = CW(keep, :); cn = cn(keep);
    if isempty(CW), error('no contrast in dsgn could be mapped onto the conditions'); end

    fprintf('\ncontrasts derived from dsgn (verify before reading the results):\n');
    fprintf('  %-32s %s\n', 'contrast', strjoin(conds, '  '));
    for c = 1:size(CW,1)
        fprintf('  %-32s %s\n', cn{c}, num2str(CW(c,:), '%+g  '));
    end
end

function e = local_entities(name)
% BIDS key-value entities from a filename, as a struct
    e = struct();
    parts = strsplit(erase(name, '.tsv'), '_');
    for i = 1:numel(parts)
        kv = strsplit(parts{i}, '-');
        if numel(kv) == 2, e.(matlab.lang.makeValidName(kv{1})) = kv{2}; end
    end
end

function [evf, how] = local_find_events(BIDSdir, sub, stem)
% Resolve the events file for one functional run.
%
% Three routes, most trustworthy first:
%
%   exact       sub-XX/func/<stem>_events.tsv
%   inherited   the BIDS inheritance principle - a file higher up the tree
%               applies to everything below it, provided every entity in its
%               name matches. A long block design with fixed timing is often
%               stored once as BIDSdir/task-<label>_events.tsv rather than
%               copied per subject and run. The most specific applicable
%               candidate wins.
%   run entity  a single events file in the subject's func dir carrying the
%               right run, used when the task label is spelled differently in
%               BIDS and in derivatives (e.g. task-emosex_movies vs
%               task-emosex), which no entity match would survive
%
    evf = ''; how = '';

    exact = fullfile(BIDSdir, sub, 'func', [stem '_events.tsv']);
    if isfile(exact), evf = exact; how = 'exact'; return, end

    want = local_entities(stem);

    % search from the most specific level of the hierarchy upwards
    levels = { fullfile(BIDSdir, sub, 'func'), fullfile(BIDSdir, sub), BIDSdir };
    if isfield(want, 'ses')
        levels = [ { fullfile(BIDSdir, sub, ['ses-' want.ses], 'func'), ...
                     fullfile(BIDSdir, sub, ['ses-' want.ses]) }, levels ];
    end

    best = ''; bestn = -1;
    for L = 1:numel(levels)
        if ~isfolder(levels{L}), continue, end
        cand = dir(fullfile(levels{L}, '*_events.tsv'));
        cand = cand(~[cand.isdir] & ~contains({cand.name}, 'noninterest'));
        for c = 1:numel(cand)
            have = local_entities(cand(c).name);
            keys = fieldnames(have);
            ok = true;
            for k = 1:numel(keys)
                if ~isfield(want, keys{k}) || ~strcmp(want.(keys{k}), have.(keys{k}))
                    ok = false; break
                end
            end
            % a candidate with more matching entities is more specific
            if ok && numel(keys) > bestn
                best = fullfile(cand(c).folder, cand(c).name); bestn = numel(keys);
            end
        end
        if ~isempty(best), evf = best; how = 'inheritance'; return, end
    end

    % last resort: right run, whatever the task label
    runtok = regexp(stem, 'run-[0-9]+', 'match', 'once');
    cand = dir(fullfile(BIDSdir, sub, 'func', '*_events.tsv'));
    cand = cand(~contains({cand.name}, 'noninterest'));
    if ~isempty(runtok), cand = cand(contains({cand.name}, runtok)); end
    if isscalar(cand)
        evf = fullfile(cand(1).folder, cand(1).name); how = 'run entity only';
    end
end

function K = local_dct(n, TR, cutoff)
    nc = fix(2*(n*TR)/cutoff + 1);
    K = ones(n,1)/sqrt(n);
    for j = 2:nc
        K(:,j) = sqrt(2/n) * cos(pi*(2*(0:n-1)'+1)*(j-1)/(2*n));
    end
end

function X = local_convolve(ev, conds, n, TR)
% SPM canonical HRF at microtime resolution 16, sampled back to volume onsets
    dt = TR/16;
    t  = (0:dt:32)';
    h  = (t.^5 .* exp(-t)/gamma(6)) - (1/6)*(t.^15 .* exp(-t)/(0.9^10*gamma(16)));
    h  = h / sum(h);
    nhi = ceil(n*TR/dt) + numel(h);
    X = zeros(n, numel(conds));
    tt = string(ev.trial_type);
    for c = 1:numel(conds)
        hi = zeros(nhi,1);
        for e = find(tt == string(conds{c}))'
            a = max(1, round(ev.onset(e)/dt) + 1);
            b = min(nhi, a + max(1, round(ev.duration(e)/dt)) - 1);
            hi(a:b) = 1;
        end
        cv = conv(hi, h);
        X(:,c) = cv(round((0:n-1)*TR/dt) + 9);
    end
end
