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
%   **'keepbetas':**
%        keep the beta images of both variants. Default false: two full
%        estimations per subject is a lot of NIfTIs to leave inside a datalad
%        subdataset, and nothing here needs them once the contrasts exist
%
%   **'keepunzipped':**
%        keep any functional images this function had to unzip. Default false,
%        so the derivatives subdataset is left exactly as it was found
%
% :Outputs:
%
%   **OUT.table:**
%        one row per subject and contrast:
%
%        - r_A_vs_original, maxdiff_A_vs_original_SE - the validation. Variant
%          A rebuilds the original fit, so it should reproduce the con images
%          already in the model directory, r > 0.999. If it does not, the
%          rebuild is not faithful and the remaining columns mean nothing; the
%          function says so rather than leaving you to notice
%        - r_A_vs_B, median_absdiff_SE, p95_absdiff_SE - the comparison itself,
%          the differences expressed in units of the contrast standard error so
%          they can be read against the leakage bound from
%          LaBGAScore_firstlevel_task_motion_diagnostics
%        - frac_voxels_flipping_p001 - proportion of in-mask voxels crossing
%          p < .001 under one variant but not the other
%
%   **OUT.outdir:**
%        where everything was written:
%
%        <outdir>/<subject>/A/                  the model refitted as it stands
%        <outdir>/<subject>/B/                  the same, corrected motion block
%        <outdir>/<subject>/diff_con_XXXX.nii   A minus B, numbered as in the
%                                               model directory
%
% :Notes:
%
% The working directory does not matter: modeldir and everything derived from
% it are absolute, and the functional images are located through the paths
% stored in SPM.xY.P. Running from the superdataset root is convenient only
% because that is where you would have run s1 to get DSGN into the workspace.
%
% The datalad datasets are left as they were found. Three things could dirty
% them, and each is handled:
%
%   - the output tree, which sits inside the first-level subdataset by default.
%     A .gitignore of '*' is written into it, so its whole contents including
%     that file are invisible to git and datalad status stays clean. Delete the
%     directory when you are done with it; nothing in it is precious
%   - the functional images. s2_fit_model gunzips them into the run
%     directories and LaBGAScore_clean_gzip_all_nii re-compresses them
%     afterwards, so SPM.mat routinely points at .nii files that no longer
%     exist. Whatever has to be unzipped to refit is removed again on the way
%     out, including if the function errors or is interrupted. Pass
%     'keepunzipped' to retain them
%   - SPM's graphics postscript, which it writes into the current directory.
%     The function works from outdir and restores the working directory when
%     it exits, so nothing is dropped into the superdataset root
%
% If the annexed content of an image is not present, the subject is skipped
% with a message naming the datalad get to run, rather than a gunzip error.
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
% LaBGAScore_firstlevel_refit_motion_comparison.m         v1.2
%
% last modified: 2026/09/03


%% PARSE OPTIONS
% -------------------------------------------------------------------------

p = inputParser;
p.addParameter('subjects', {}, @iscell);
p.addParameter('outdir', '', @(x) ischar(x) || isstring(x));
p.addParameter('nmotion', 24, @(x) isnumeric(x) && ismember(x, [12 24]));
p.addParameter('contrasts', [], @isnumeric);
p.addParameter('mask', '', @(x) ischar(x) || isstring(x));
p.addParameter('keepbetas', false, @islogical);
p.addParameter('keepunzipped', false, @islogical);
p.parse(varargin{:});
opt = p.Results;

if isempty(which('spm')), error('SPM12 is not on the Matlab path'); end

modeldir = char(modeldir);
if ~isfolder(modeldir), error('modeldir does not exist: %s', modeldir); end
dm = dir(modeldir); modeldir = dm(1).folder;      % absolute, since we cd below

outdir = char(opt.outdir);
if isempty(outdir), outdir = fullfile(modeldir, 'motion_refit_comparison'); end
if ~isfolder(outdir), mkdir(outdir); end
do = dir(outdir); outdir = do(1).folder;

% Keep the dataset clean. The output tree is bulky and regenerable, so it is
% made invisible to git rather than left for datalad status to report. A
% .gitignore of '*' ignores the directory's whole contents, itself included,
% which works whether or not outdir sits inside a subdataset.
gi = fullfile(outdir, '.gitignore');
if ~isfile(gi)
    fid = fopen(gi, 'w');
    if fid > 0
        fprintf(fid, ['# written by LaBGAScore_firstlevel_refit_motion_comparison\n' ...
                      '# everything here is regenerable, and is kept out of git so the\n' ...
                      '# subdataset stays clean\n*\n']);
        fclose(fid);
    end
end

% SPM writes its graphics postscript into the CURRENT directory, which would
% otherwise land in the superdataset root. Work from outdir instead, and go
% back afterwards however this function exits.
origwd = pwd;
cleanupWd = onCleanup(@() cd(origwd));
cd(outdir);

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

    % ---- make sure the functional images are actually on disk -------------
    % s2_fit_model gunzips the smoothed images into the run directories, but
    % LaBGAScore_clean_gzip_all_nii re-compresses them afterwards to keep the
    % dataset small, so SPM.xY.P routinely points at .nii files that no longer
    % exist. Unzip what is needed and put it back afterwards.
    try
        unz = local_ensure_images(S);
    catch ME
        warning('%s: %s', subs{s}, ME.message); continue
    end
    cleanupImgs = onCleanup(@() local_remove_files(unz, opt.keepunzipped));

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
        if ~opt.keepbetas
            % two full estimations per subject is a lot of NIfTIs to leave
            % inside a datalad subdataset, and the betas are not needed once the
            % contrasts exist; ResMS and mask are, so they stay
            delete(fullfile(vd, 'beta_*.nii'));
        end
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

        % contrasts are renumbered 1..numel(ci) inside the variant directories,
        % while the original model numbers them by their index in SPM.xCon
        cA = spm_read_vols(spm_vol(fullfile(vdirs.A, sprintf('con_%04d.nii', k))));
        cB = spm_read_vols(spm_vol(fullfile(vdirs.B, sprintf('con_%04d.nii', k))));
        tA = spm_read_vols(spm_vol(fullfile(vdirs.A, sprintf('spmT_%04d.nii', k))));
        tB = spm_read_vols(spm_vol(fullfile(vdirs.B, sprintf('spmT_%04d.nii', k))));

        a = cA(:); b = cB(:); ta = tA(:); tb = tB(:);
        keep = m & isfinite(a) & isfinite(b) & isfinite(ta) & isfinite(tb);

        % difference in units of the contrast SE, taken from variant A
        cc = S.xCon(ci(k)).c;
        seA = sqrt(resA(:) * (cc' * S.xX.Bcov * cc));
        dse = abs(a - b) ./ max(seA, eps);

        dfe = S.xX.erdf;
        thr = spm_invTcdf(1 - 0.001, dfe);
        flip = xor(abs(ta) > thr, abs(tb) > thr);

        % validation: variant A rebuilds the original fit, so its con image
        % should be all but identical to the one already in the model directory.
        % Anything below ~0.999 means the rebuild is not reproducing the
        % original model and the A-versus-B comparison should not be read.
        rOrig = NaN; dOrig = NaN;
        origfile = fullfile(modeldir, subs{s}, sprintf('con_%04d.nii', ci(k)));
        if isfile(origfile)
            co = spm_read_vols(spm_vol(origfile)); o = co(:);
            ko = keep & isfinite(o);
            if any(ko)
                rOrig = corr(a(ko), o(ko));
                dOrig = max(abs(a(ko) - o(ko)) ./ max(seA(ko), eps));
            end
        end

        rows(end+1,:) = { subs{s}, S.xCon(ci(k)).name, rOrig, dOrig, corr(a(keep), b(keep)), ...
            median(dse(keep)), prctile(dse(keep), 95), mean(flip(keep)) }; %#ok<AGROW>

        % write the difference image for inspection, named by the ORIGINAL
        % contrast number so it lines up with the model directory
        V = spm_vol(fullfile(vdirs.A, sprintf('con_%04d.nii', k)));
        V.fname = fullfile(outdir, subs{s}, sprintf('diff_con_%04d.nii', ci(k)));
        V.descrip = sprintf('A minus B: %s', S.xCon(ci(k)).name);
        spm_write_vol(V, cA - cB);
    end
end

OUT.table = cell2table(rows, 'VariableNames', ...
    {'subject','contrast','r_A_vs_original','maxdiff_A_vs_original_SE', ...
     'r_A_vs_B','median_absdiff_SE','p95_absdiff_SE','frac_voxels_flipping_p001'});
OUT.outdir = outdir;

fprintf('\n');
disp(OUT.table);

bad = OUT.table.r_A_vs_original < 0.999 & ~isnan(OUT.table.r_A_vs_original);
if any(bad)
    warning(['VALIDATION FAILED for %d of %d rows: variant A does not reproduce the ' ...
             'original con images (lowest r = %.4f). The rebuild is not reproducing the ' ...
             'original model, so the A-versus-B columns should NOT be interpreted. ' ...
             'Check that the images SPM.xY.P points at are the ones the model was fitted on.'], ...
             sum(bad), height(OUT.table), min(OUT.table.r_A_vs_original(bad)));
elseif all(isnan(OUT.table.r_A_vs_original))
    warning('no original con images found in %s, so the rebuild could not be validated', modeldir);
else
    fprintf('validation passed: variant A reproduces the original con images (min r = %.5f)\n', ...
        min(OUT.table.r_A_vs_original));
end

fprintf(['\nr_A_vs_B near 1 with median_absdiff_SE near 0 means the parameterisation does\n' ...
         'not matter here. frac_voxels_flipping_p001 is the proportion of in-mask voxels\n' ...
         'crossing p < .001 under one variant but not the other.\n' ...
         '\nWritten to %s\n' ...
         '  <subject>/A/          the model refitted as it stands\n' ...
         '  <subject>/B/          the same model with the corrected motion block\n' ...
         '  <subject>/diff_con_XXXX.nii   A minus B, numbered as in the model directory\n'], outdir);
if ~opt.keepbetas
    fprintf('\nbeta images were deleted after contrast estimation; pass ''keepbetas'', true to keep them\n');
end
fprintf(['\nThe datasets are as they were found: the output tree carries a .gitignore of\n' ...
         '''*'' so datalad status stays clean, any functional images unzipped to refit have\n' ...
         'been removed, and the working directory is back at %s. Verify with:\n' ...
         '  datalad status -r\n'], origwd);

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

function unz = local_ensure_images(S)
% Every volume in SPM.xY.P must exist before SPM can map the files. Where the
% .nii has been re-compressed since the model was fitted, unzip it and return
% the list so it can be removed again afterwards.
    unz = {};
    files = unique(regexprep(cellstr(S.xY.P), ',\d+$', ''), 'stable');
    missing = {};
    for k = 1:numel(files)
        f = files{k};
        if isfile(f), continue, end
        gz = [f '.gz'];
        if ~isfile(gz), missing{end+1} = f; continue, end %#ok<AGROW>
        d = dir(gz);
        if isempty(d) || d.bytes == 0
            % a git-annex pointer whose content has not been fetched
            error(['%s is present but empty, which usually means the annexed ' ...
                   'content is not here. Run: datalad get %s'], gz, gz);
        end
        gunzip(gz, fileparts(f));
        if ~isfile(f)
            error('gunzip of %s did not produce %s', gz, f);
        end
        unz{end+1} = f; %#ok<AGROW>
    end
    if ~isempty(missing)
        error(['%d functional file(s) referenced by SPM.mat are missing, with no .gz ' ...
               'alongside, the first being %s. The model cannot be refitted without ' ...
               'the data it was fitted on.'], numel(missing), missing{1});
    end
    if ~isempty(unz)
        fprintf('  unzipped %d functional file(s) that had been re-compressed\n', numel(unz));
    end
end

function local_remove_files(files, keep)
    if keep || isempty(files), return, end
    for k = 1:numel(files)
        if isfile(files{k}), delete(files{k}); end
    end
    fprintf('  removed %d unzipped functional file(s), leaving the dataset as it was\n', numel(files));
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
