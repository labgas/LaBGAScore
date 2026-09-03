function OUT = LaBGAScore_firstlevel_refit_motion_comparison(modeldir, varargin)
% Refit first-level models with corrected motion regressors and compare con images
%
% :Usage:
% ::
%
%     OUT = LaBGAScore_firstlevel_refit_motion_comparison(modeldir, 'subjects', {'sub-001','sub-003'})
%
% Measures what the motion-regressor parameterisation actually does to a
% study's contrast estimates, which no amount of design-level reasoning can
% settle. For each subject it fits the SAME model twice, changing one thing:
%
%   variant A   the motion block exactly as it sits in the existing SPM.mat
%   variant B   the motion block rebuilt as z-scored [parameters, first
%               derivatives] (+ quadratics), i.e. the Friston-style expansion
%               the script headers describe
%
% Everything else - the task regressors, the high-pass filter, the AR model,
% the spike and CSF regressors, the masking threshold, the contrast vectors -
% is taken from the existing SPM.mat and held identical, so any difference in
% the con images is attributable to the motion columns alone.
%
% Nothing is written into the existing model directory. Both variants are
% estimated in a scratch tree under 'outdir'.
%
% :Inputs:
%
%   **modeldir:**
%        first-level model directory holding sub-*/SPM.mat, i.e. DSGN.modeldir
%
% :Optional Inputs:
%
%   **'subjects':**
%        cell array of subject directory names. Default: the first three
%        found, which is normally enough to see whether this matters
%
%   **'outdir':**
%        where to estimate the two variants. Default:
%        <modeldir>/motion_refit_comparison
%
%   **'nmotion':**
%        number of leading motion columns in the multiple-regressor block:
%        24 with quadratics, 12 without. Default 24. The value is verified
%        against the design before anything is estimated (see Notes), so a
%        wrong setting errors rather than producing a misleading comparison
%
%   **'contrasts':**
%        indices into SPM.xCon to compare. Default: all
%
%   **'mask':**
%        image to restrict the comparison to. Default: the model's own mask.nii
%
% :Outputs:
%
%   **OUT:**
%        struct with a per-subject-per-contrast table of the correlation
%        between the two variants' con images, the median absolute difference
%        expressed in units of the contrast's standard error, and the
%        proportion of in-mask voxels whose t crosses p < .001 under one
%        variant but not the other. Difference images are written alongside
%        the estimates for inspection
%
% :Notes:
%
% BEFORE estimating anything, the function rebuilds what it believes the
% fitted motion block to be from the run's own confound file and checks it
% against the first 'nmotion' columns actually present in SPM.Sess(i).C.C. If
% they do not match to near machine precision the run is refused, because a
% mismatch means either 'nmotion' is wrong or the regressors were not built by
% the LaBGAScore pipeline, and in both cases variant B would not be the
% controlled comparison it claims to be.
%
% Estimating two variants for one subject costs roughly twice a normal
% first-level fit. Three subjects is a sensible first pass.
%
% :See also: LaBGAScore_firstlevel_task_motion_diagnostics,
% LaBGAScore_firstlevel_s2_fit_model
%
% -------------------------------------------------------------------------
%
% authors: Lukas Van Oudenhove
%
% date: KU Leuven, September, 2026
%
% -------------------------------------------------------------------------
%
% LaBGAScore_firstlevel_refit_motion_comparison.m         v1.0
%
% last modified: 2026/09/02


%% PARSE OPTIONS
% -------------------------------------------------------------------------

p = inputParser;
p.addParameter('subjects', {}, @iscell);
p.addParameter('outdir', '', @(x) ischar(x) || isstring(x));
p.addParameter('nmotion', 24, @(x) isnumeric(x) && ismember(x, [12 24]));
p.addParameter('contrasts', [], @isnumeric);
p.addParameter('mask', '', @(x) ischar(x) || isstring(x));
p.parse(varargin{:});
opt = p.Results;

if isempty(which('spm')), error('SPM12 is not on the Matlab path'); end

modeldir = char(modeldir);
outdir = char(opt.outdir);
if isempty(outdir), outdir = fullfile(modeldir, 'motion_refit_comparison'); end
if ~isfolder(outdir), mkdir(outdir); end

d = dir(fullfile(modeldir, 'sub-*'));
subs = {d([d.isdir]).name}';
if ~isempty(opt.subjects), subs = intersect(subs, opt.subjects, 'stable'); end
if isempty(subs), error('no subject directories found in %s', modeldir); end
if isempty(opt.subjects), subs = subs(1:min(3, numel(subs))); end

spm('defaults', 'FMRI');
spm_jobman('initcfg');

rows = {};


%% LOOP OVER SUBJECTS
% -------------------------------------------------------------------------

for s = 1:numel(subs)

    fprintf('\n=== %s ===\n', subs{s});

    spmfile = fullfile(modeldir, subs{s}, 'SPM.mat');
    if ~isfile(spmfile), warning('no SPM.mat for %s, skipping', subs{s}); continue, end
    L = load(spmfile); S = L.SPM;

    % ---- rebuild the two motion blocks, and verify the first one ----------
    nsess = numel(S.Sess);
    regsA = cell(nsess,1); regsB = cell(nsess,1); ok = true;

    for i = 1:nsess
        C = S.Sess(i).C.C;                       % the multiple-regressor block as fitted
        if size(C,2) < opt.nmotion
            warning('%s session %d has %d regressors, fewer than nmotion=%d, skipping subject', ...
                subs{s}, i, size(C,2), opt.nmotion); ok = false; break
        end

        R = local_raw_motion(S, i);
        if isempty(R), ok = false; break, end

        [fitA, fitB] = local_motion_blocks(R, opt.nmotion);

        % refuse to proceed unless the reconstruction matches what was fitted
        r = zeros(1, opt.nmotion);
        for k = 1:opt.nmotion
            a = C(:,k); b = fitA(:,k);
            if std(a) < eps || std(b) < eps, r(k) = double(std(a) < eps && std(b) < eps);
            else, r(k) = abs(corr(a, b)); end
        end
        if min(r) < 0.999
            warning(['%s session %d: rebuilt motion block does not match the design ' ...
                     '(worst |r| = %.4f). Check nmotion, or whether these regressors ' ...
                     'came from the LaBGAScore pipeline. Skipping subject.'], subs{s}, i, min(r));
            ok = false; break
        end

        regsA{i} = [fitA C(:, opt.nmotion+1:end)];   % as fitted
        regsB{i} = [fitB C(:, opt.nmotion+1:end)];   % intended expansion, same spikes and CSF
    end
    if ~ok, continue, end

    % ---- estimate both variants ------------------------------------------
    vdirs = struct();
    for v = {'A','B'}
        tag = v{1};
        vd = fullfile(outdir, subs{s}, tag);
        if isfolder(vd), rmdir(vd, 's'); end
        mkdir(vd);
        regs = regsA; if strcmp(tag,'B'), regs = regsB; end

        regfiles = cell(nsess,1);
        for i = 1:nsess
            R = regs{i};
            regfiles{i} = fullfile(vd, sprintf('noise_regs_sess%02d.mat', i));
            save(regfiles{i}, 'R');
        end

        local_specify_estimate(S, vd, regfiles);
        local_contrasts(S, vd, opt.contrasts);
        vdirs.(tag) = vd;
        fprintf('  variant %s estimated\n', tag);
    end

    % ---- compare ---------------------------------------------------------
    maskfile = char(opt.mask);
    if isempty(maskfile), maskfile = fullfile(vdirs.A, 'mask.nii'); end
    mv = spm_read_vols(spm_vol(maskfile)); m = mv(:) > 0;

    ci = opt.contrasts; if isempty(ci), ci = 1:numel(S.xCon); end
    resA = spm_read_vols(spm_vol(fullfile(vdirs.A, 'ResMS.nii')));

    for k = 1:numel(ci)
        cA = spm_read_vols(spm_vol(fullfile(vdirs.A, sprintf('con_%04d.nii', ci(k)))));
        cB = spm_read_vols(spm_vol(fullfile(vdirs.B, sprintf('con_%04d.nii', ci(k)))));
        tA = spm_read_vols(spm_vol(fullfile(vdirs.A, sprintf('spmT_%04d.nii', ci(k)))));
        tB = spm_read_vols(spm_vol(fullfile(vdirs.B, sprintf('spmT_%04d.nii', ci(k)))));

        a = cA(:); b = cB(:); ta = tA(:); tb = tB(:);
        keep = m & isfinite(a) & isfinite(b) & isfinite(ta) & isfinite(tb);

        % difference in units of the contrast SE, taken from variant A
        cc = S.xCon(ci(k)).c;
        seA = sqrt(resA(:) * (cc' * S.xX.Bcov * cc));
        dse = abs(a - b) ./ max(seA, eps);

        dfe = S.xX.erdf;
        thr = spm_invTcdf(1 - 0.001, dfe);
        flip = xor(abs(ta) > thr, abs(tb) > thr);

        rows(end+1,:) = { subs{s}, S.xCon(ci(k)).name, corr(a(keep), b(keep)), ...
            median(dse(keep)), prctile(dse(keep), 95), mean(flip(keep)) }; %#ok<AGROW>

        % write the difference image for inspection
        V = spm_vol(fullfile(vdirs.A, sprintf('con_%04d.nii', ci(k))));
        V.fname = fullfile(outdir, subs{s}, sprintf('diff_con_%04d.nii', ci(k)));
        V.descrip = 'variant A minus variant B';
        spm_write_vol(V, cA - cB);
    end
end

OUT.table = cell2table(rows, 'VariableNames', ...
    {'subject','contrast','r_con','median_absdiff_SE','p95_absdiff_SE','frac_voxels_flipping_p001'});
OUT.outdir = outdir;

fprintf('\n');
disp(OUT.table);
fprintf(['\nr_con near 1 and median_absdiff_SE near 0 means the parameterisation does not\n' ...
         'matter here. frac_voxels_flipping_p001 is the proportion of in-mask voxels that\n' ...
         'cross p < .001 under one variant but not the other. Difference images are in %s\n'], outdir);

end % main


%% LOCAL FUNCTIONS
% -------------------------------------------------------------------------

function R = local_raw_motion(S, sess)
% the six realignment parameters for one session, from the confound file that
% LaBGAScore_firstlevel_s2_fit_model copied next to the functional images
    R = [];
    rows = local_session_rows(S, sess);
    img = strtrim(S.xY.P(rows(1), :));
    rundir = fileparts(regexprep(img, ',\d+$', ''));
    c = dir(fullfile(rundir, '*desc-confounds_timeseries.tsv'));
    if isempty(c), c = dir(fullfile(fileparts(rundir), '*desc-confounds_timeseries.tsv')); end
    if isempty(c)
        warning('no confound file found near %s', rundir); return
    end
    T = readtable(fullfile(c(1).folder, c(1).name), 'FileType','text', 'Delimiter','\t', ...
        'TreatAsEmpty','n/a', 'VariableNamingRule','preserve');
    mot = {'trans_x','trans_y','trans_z','rot_x','rot_y','rot_z'};
    if ~all(ismember(mot, T.Properties.VariableNames))
        warning('%s lacks the six realignment parameters', c(1).name); R = []; return
    end
    R = T{:, mot}; R(~isfinite(R)) = 0;
    if size(R,1) ~= numel(rows)
        warning('confound file has %d rows but session %d has %d scans', size(R,1), sess, numel(rows));
        R = []; return
    end
end

function rows = local_session_rows(S, sess)
    ns = S.nscan; st = cumsum([0 ns(:)']);
    rows = (st(sess)+1) : st(sess+1);
end

function [fitA, fitB] = local_motion_blocks(R, nmotion)
% fitA reproduces what the pipeline fitted (first and second derivatives),
% fitB is the intended expansion (parameters and first derivatives, z-scored)
    g = @(X) cell2mat(arrayfun(@(c) gradient(X(:,c)), 1:size(X,2), 'UniformOutput', false));
    zs = @(X) (X - mean(X)) ./ max(std(X), eps);
    d1 = g(R); d2 = g(d1);
    fitA = [d1 d2];
    fitB = zs([R d1]);
    if nmotion == 24
        fitA = [fitA fitA.^2];
        fitB = [fitB fitB.^2];
    end
end

function local_specify_estimate(S, vd, regfiles)
% rebuild the model exactly as specified in S, but pointing at new regressors
    clear matlabbatch
    matlabbatch{1}.spm.stats.fmri_spec.dir = {vd};
    matlabbatch{1}.spm.stats.fmri_spec.timing.units   = S.xBF.UNITS;
    matlabbatch{1}.spm.stats.fmri_spec.timing.RT      = S.xY.RT;
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t  = S.xBF.T;
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = S.xBF.T0;

    for i = 1:numel(S.Sess)
        rows = local_session_rows(S, i);
        matlabbatch{1}.spm.stats.fmri_spec.sess(i).scans = cellstr(S.xY.P(rows, :));
        for j = 1:numel(S.Sess(i).U)
            U = S.Sess(i).U(j);
            matlabbatch{1}.spm.stats.fmri_spec.sess(i).cond(j).name     = U.name{1};
            matlabbatch{1}.spm.stats.fmri_spec.sess(i).cond(j).onset    = U.ons;
            matlabbatch{1}.spm.stats.fmri_spec.sess(i).cond(j).duration = U.dur;
            matlabbatch{1}.spm.stats.fmri_spec.sess(i).cond(j).tmod     = 0;
            matlabbatch{1}.spm.stats.fmri_spec.sess(i).cond(j).orth     = 1;
            pm = struct('name', {}, 'param', {}, 'poly', {});
            if isfield(U, 'P') && ~isempty(U.P) && ~strcmp(U.P(1).name, 'none')
                for q = 1:numel(U.P)
                    pm(q).name  = U.P(q).name;
                    pm(q).param = U.P(q).P;
                    pm(q).poly  = U.P(q).h;
                end
            end
            matlabbatch{1}.spm.stats.fmri_spec.sess(i).cond(j).pmod = pm;
        end
        matlabbatch{1}.spm.stats.fmri_spec.sess(i).multi     = {''};
        matlabbatch{1}.spm.stats.fmri_spec.sess(i).regress   = struct('name', {}, 'val', {});
        matlabbatch{1}.spm.stats.fmri_spec.sess(i).multi_reg = regfiles(i);
        matlabbatch{1}.spm.stats.fmri_spec.sess(i).hpf       = S.xX.K(i).HParam;
    end

    matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
    derivs = [0 0];
    if isfield(S.xBF, 'order') && S.xBF.order == 3, derivs = [1 1];
    elseif isfield(S.xBF, 'order') && S.xBF.order == 2, derivs = [1 0]; end
    matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = derivs;
    matlabbatch{1}.spm.stats.fmri_spec.volt    = S.xBF.Volterra;
    matlabbatch{1}.spm.stats.fmri_spec.global  = 'None';
    matlabbatch{1}.spm.stats.fmri_spec.mthresh = S.xM.gMT;
    matlabbatch{1}.spm.stats.fmri_spec.mask    = {''};
    matlabbatch{1}.spm.stats.fmri_spec.cvi     = S.xVi.form;

    matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = ...
        cfg_dep('fMRI model specification: SPM.mat File', ...
        substruct('.','val','{}',{1},'.','val','{}',{1},'.','val','{}',{1}), substruct('.','spmmat'));
    matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
    matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

    spm_jobman('run', matlabbatch);
end

function local_contrasts(S, vd, ci)
% reuse the original contrast vectors: the design has the same columns in the
% same order, only the motion values differ
    if isempty(ci), ci = 1:numel(S.xCon); end
    clear matlabbatch
    matlabbatch{1}.spm.stats.con.spmmat = {fullfile(vd, 'SPM.mat')};
    for k = 1:numel(ci)
        matlabbatch{1}.spm.stats.con.consess{k}.tcon.name    = S.xCon(ci(k)).name;
        matlabbatch{1}.spm.stats.con.consess{k}.tcon.weights = S.xCon(ci(k)).c';
        matlabbatch{1}.spm.stats.con.consess{k}.tcon.sessrep = 'none';
    end
    matlabbatch{1}.spm.stats.con.delete = 1;
    spm_jobman('run', matlabbatch);
end
