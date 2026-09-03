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
%        BIDS directory holding sub-*/func/*_events.tsv
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
%   **'contrasts':**
%        [nContrast x nCondition] weight matrix. Default: one contrast per
%        condition against implicit baseline
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
% LaBGAScore_firstlevel_task_motion_diagnostics.m         v1.0
%
% last modified: 2026/09/02


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
if isempty(conds)
    seen = string([]);
    for s = 1:numel(subs)
        e = dir(fullfile(BIDSdir, subs{s}, 'func', [taskpat '_events.tsv']));
        for f = 1:numel(e)
            ev = readtable(fullfile(e(f).folder, e(f).name), 'FileType','text', ...
                'Delimiter','\t', 'VariableNamingRule','preserve');
            seen = union(seen, unique(string(ev.trial_type)));
        end
    end
    conds = cellstr(seen);
end
nC = numel(conds);

CW = opt.contrasts; cn = opt.contrastnames;
if isempty(CW), CW = eye(nC); end
if isempty(cn), cn = arrayfun(@(k) sprintf('contrast %d', k), 1:size(CW,1), 'UniformOutput', false); end


%% LOOP OVER SUBJECTS AND RUNS
% -------------------------------------------------------------------------

rows = {}; Lk = []; Lsub = {}; Lidx = [];
mot = {'trans_x','trans_y','trans_z','rot_x','rot_y','rot_z'};

for s = 1:numel(subs)

    cfs = dir(fullfile(derivdir, subs{s}, 'func', [taskpat '_desc-confounds_timeseries.tsv']));

    for f = 1:numel(cfs)

        % match the events file to this run by stripping the fMRIPrep suffix
        stem = erase(cfs(f).name, '_desc-confounds_timeseries.tsv');
        evf  = fullfile(BIDSdir, subs{s}, 'func', [stem '_events.tsv']);
        if ~isfile(evf)
            % the task label is not always spelled the same way in BIDS and in
            % derivatives (e.g. task-emosex_movies vs task-emosex), so fall
            % back to matching on subject and run entity alone
            runtok = regexp(stem, 'run-[0-9]+', 'match', 'once');
            cand = dir(fullfile(BIDSdir, subs{s}, 'func', '*_events.tsv'));
            cand = cand(~contains({cand.name}, 'noninterest'));
            if ~isempty(runtok)
                cand = cand(contains({cand.name}, runtok));
            end
            if numel(cand) == 1
                evf = fullfile(cand(1).folder, cand(1).name);
            elseif isempty(cand)
                warning('no events file for %s, skipping', stem); continue
            else
                warning('%d candidate events files for %s, skipping', numel(cand), stem); continue
            end
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

        Kr = @(v) v - K*(K\v);
        for c = 1:nC
            x = X(:,c) - mean(X(:,c));
            if std(x) < eps, continue; end
            xr = x - QN*(QN'*x);            % task variance surviving the fitted nuisance model
            rows(end+1,:) = { subs{s}, stem, conds{c}, ...
                norm(QM'*xr)^2 / max(norm(xr)^2, eps), ...
                corr(Kr(x), Kr(fd)), num2cell(corr(Kr(x), Kr(R)))  }; %#ok<AGROW>
        end

        % per-contrast leakage bound
        for k = 1:size(CW,1)
            w = zeros(1, size(Df,2)); w(1:nC) = CW(k,:);
            r = w * Dpf;
            Lk(end+1,1) = norm(r*QA) / norm(r); %#ok<AGROW>
            Lsub{end+1,1} = subs{s}; Lidx(end+1,1) = k; %#ok<AGROW>
        end
    end
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

end % main


%% LOCAL FUNCTIONS
% -------------------------------------------------------------------------

function d = local_grad(X)
    d = cell2mat(arrayfun(@(c) gradient(X(:,c)), 1:size(X,2), 'UniformOutput', false));
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
