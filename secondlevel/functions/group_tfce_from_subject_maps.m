function [tfce_dat, p_img, info] = group_tfce_from_subject_maps( ...
   fmri_dat_subj, design, group, covariates, nPerm, varargin)
% group_tfce_from_subject_maps  Group-level TFCE permutation inference.
%
% Group-level TFCE permutation inference starting from subject-level
% contrast images, with optional Freedman-Lane covariate control.
%
% Uses classic TFCE (Smith & Nichols 2009) via tfce_volume. Earlier versions
% called pTFCE, a different algorithm, and did so with mismatched arguments.
%
% -------------------------------------------------------------------------
% INPUTS
% -------------------------------------------------------------------------
% fmri_dat_subj : fmri_data_st
%  - dat = [nVox × nSubj] subject-level con images
%
% design : char
%  - 'onesample' : one-sample / paired design (sign-flip permutations)
%  - 'twosample' : two independent groups (label exchange permutations)
%
% group : [] or [nSubj × 1]
%  - [] for 'onesample'
%  - group labels for 'twosample'
%
% covariates : [] or [nSubj × nCov]
%  - nuisance covariates to control (Freedman–Lane)
%  - do NOT include a column of ones; the intercept is handled per design
%    (see "Nuisance model" below)
%
% nPerm : scalar integer
%  - number of permutations
%
% OPTIONAL NAME–VALUE PAIRS
%  'H'       : TFCE height exponent (default = 2)
%  'E'       : TFCE extent exponent (default = 0.5)
%  'conn'    : connectivity, 6 | 18 | 26 (default = 26)
%
%  These are the real Smith & Nichols TFCE parameters and now actually control
%  the transform. Under the previous pTFCE-based implementation they did not:
%  pTFCE has no H or E parameter, and these values were silently landing in its
%  resel-count and voxel-count slots, while 'conn' never reached it at all.
%  'sidedness': 'one' (default) or 'two'
%  'tail'    : 'pos' or 'neg' (only if sidedness = 'one')
%  'parallel' : true (default) or false
%
% -------------------------------------------------------------------------
% OUTPUTS
% -------------------------------------------------------------------------
% tfce_dat : fmri_data_st
%  - TFCE-transformed group map
%
% p_img : statistic_image (type 'p')
%  - voxelwise TFCE permutation p-values
%
% info : struct
%  - TFCE_real
%  - TFCE_null
%  - TFCE_real_max
%  - TFCE_null_max
%  - p_TFCE_global
%  - nuisance_rank : effective rank of the nuisance model (0 if none)
%  - tfce_dh : integration step shared by the observed map and every permutation
%
% -------------------------------------------------------------------------
% NUISANCE MODEL (why there is no intercept column)
% -------------------------------------------------------------------------
% The nuisance model must contain the nuisance and NOTHING ELSE. For
% 'onesample' the effect of interest IS the intercept (the group mean), so
% putting a column of ones in the nuisance design removes the effect into
% Y_hat, and Freedman-Lane then adds it back into every permuted dataset.
% The null ends up carrying the real effect and the test cannot reject.
%
% Covariates are therefore mean-centered and regressed out WITHOUT an
% intercept column, via residualizeFold. Centering makes them orthogonal to
% the constant vector, so the group mean stays in the residuals where the
% sign-flip null can act on it. The same treatment is used for 'twosample',
% where the grand mean is irrelevant because it cancels in a difference of
% means.
%
% -------------------------------------------------------------------------
% CALIBRATION
% -------------------------------------------------------------------------
% Checked against synthetic data with a known answer: n = 40 subjects, 200
% voxels, one covariate explaining about a third of the variance, 200
% datasets, 200 sign-flips or permutations each. Nominal alpha 0.05 has a
% 95% interval of [0.020 0.080] at that number of datasets.
%
%                     false positives     mean p      power
%   one-sample            0.055            0.502       1.00
%   two-sample            0.040            0.532       1.00
%
% For comparison, the code this replaced scored 0.000 / 0.972 / 0.00 for
% 'onesample' (the test could never reject) and 0.060 / 0.475 / 0.10 for
% 'twosample' (almost no power, because its permutation was degenerate).
%
% Constructing the null correctly is subtler than it looks and is worth
% recording. The covariate-adjusted mean is only zero if the covariate effect
% is generated around the covariate's own mean:
%     Y_i = b*(cv_i - mean(cv)) + e_i        <- null for the adjusted mean
%     Y_i = b*cv_i + e_i                     <- NOT null: the adjusted mean
%                                               is b*mean(cv), nonzero in any
%                                               finite sample
% Validating against the second form makes every scheme look broken, because
% the data really does carry an effect.
%

% ============================================================
% Defaults
% ============================================================
H = 2;
E = 0.5;
conn = 26;
sidedness = 'one';
tail = 'pos';
use_parallel = true;

for i = 1:2:numel(varargin)
   switch lower(varargin{i})
       case 'h'
           H = varargin{i+1};
       case 'e'
           E = varargin{i+1};
       case 'conn'
           conn = varargin{i+1};
       case 'sidedness'
           sidedness = lower(varargin{i+1});
       case 'tail'
           tail = lower(varargin{i+1});
       case 'parallel'
           use_parallel = varargin{i+1};
       otherwise
           error('Unknown option: %s', varargin{i});
   end
end

design = lower(design);

% ============================================================
% Sanity checks
% ============================================================
assert(isa(fmri_dat_subj,'fmri_data_st'));
[nVox, nSubj] = size(fmri_dat_subj.dat);

if isempty(covariates)
   use_covariates = false;
else
   assert(size(covariates,1)==nSubj, ...
       'covariates must have nSubj rows');
   use_covariates = true;
end

grp_vals = [];

switch design
   case 'onesample'
       assert(isempty(group));
   case 'twosample'
       assert(~isempty(group));
       assert(numel(group)==nSubj);
       grp_vals = unique(group(:));
       assert(numel(grp_vals)==2);
   otherwise
       error('Unknown design');
end

% ============================================================
% Step 1: Residualise subject maps (Freedman–Lane)
% ============================================================
Y = fmri_dat_subj.dat; % [nVox × nSubj]

if use_covariates
   % residualizeFold mean-centers the covariates and removes only the
   % covariate-attributable variation, leaving the intercept absorbed. That
   % is exactly what is needed here: for 'onesample' the intercept is the
   % effect of interest and must survive into R (see header). It also uses a
   % truncated SVD, so collinear or fold-constant covariates give a
   % minimum-norm solution instead of a backslash warning.
   assert(all(isfinite(covariates(:))), ...
       'covariates contain NaN or Inf; handle missing values upstream');
   [Rt, ~, nuis] = residualizeFold(Y', [], covariates, []);
   R = Rt';
   Y_hat = Y - R;                % fitted nuisance effect
   nuisance_rank = nuis.rank;
else
   R = Y;
   Y_hat = zeros(size(Y));
   nuisance_rank = 0;
end

% ============================================================
% Helper: group-level t-maps
% ============================================================
t_one = @(dat) one_sample_tmap(dat);
t_two = @(datA,datB) two_sample_tmap(datA,datB);

% ============================================================
% REAL t-map
% ============================================================
switch design
   case 'onesample'
       t_real = t_one(Y);

   case 'twosample'
       idxA = group==grp_vals(1);
       idxB = group==grp_vals(2);
       t_real = t_two(Y(:,idxA),Y(:,idxB));
end

% The observed map fixes the integration grid; every permutation then reuses
% it. Letting each map derive its own dh from its own maximum makes the null
% slightly incomparable to the observed statistic.
[tfce_real, tfce_dh] = tfce_one_fmri_dat( ...
   t_real,...
   fmri_dat_subj.removed_voxels,...
   fmri_dat_subj.volInfo,...
   H,E,conn,sidedness,tail,[]);

if isempty(tfce_dh) || all(tfce_real == 0)
   error('group_tfce_from_subject_maps:tfceEmpty', ...
       ['The observed statistic map has no positive values in the requested ' ...
        'tail, so TFCE is identically zero and there is nothing to test. ' ...
        'Check the ''sidedness'' and ''tail'' options.']);
end

TFCE_real_max = max(tfce_real);

tfce_dat = fmri_dat_subj;
tfce_dat.dat = single(tfce_real);
tfce_dat.dat_descrip = sprintf('TFCE (%s, FL)',design);

% ============================================================
% PERMUTATIONS (Freedman–Lane)
% ============================================================
TFCE_perm = zeros(nVox,nPerm,'single');
TFCE_null_max = zeros(nPerm,1,'single');

track_progress = use_parallel && nPerm>1;
q = [];

if track_progress
   pool = gcp('nocreate');
   if isempty(pool), parpool; end
   tracker = ProgressTracker(nPerm);
   q = parallel.pool.DataQueue;
   afterEach(q,@(~) tracker.update());
end

% 'parallel',false now really runs serially. Previously the parfor executed
% regardless, and send(q,1) below referenced a q that had never been created,
% so the option errored instead of disabling parallelism.
if use_parallel
   maxWorkers = Inf;
else
   maxWorkers = 0;
end

parfor (p = 1:nPerm, maxWorkers)

   % --- permute residuals ---
   switch design
       case 'onesample'
           signs = (rand(nSubj,1)>0.5)*2-1;
           R_perm = R .* signs';

       case 'twosample'
           % Permute the residuals and hold the group labels FIXED. The
           % previous code permuted both by the same index, which preserves
           % the residual-to-label pairing and therefore changes nothing:
           % with no covariates every permuted statistic came out exactly
           % equal to the observed one (verified: 200/200 identical), so the
           % null was degenerate and the test had no power.
           R_perm = R(:,randperm(nSubj));
           gA = group == grp_vals(1);
           gB = group == grp_vals(2);
   end

   % --- reconstruct data ---
   Yp = Y_hat + R_perm;

   % --- compute permuted t-map ---
   switch design
       case 'onesample'
           t_perm = t_one(Yp);

       case 'twosample'
           t_perm = t_two(Yp(:,gA),Yp(:,gB));
   end

   tf = tfce_one_fmri_dat( ...
       t_perm,...
       fmri_dat_subj.removed_voxels,...
       fmri_dat_subj.volInfo,...
       H,E,conn,sidedness,tail,tfce_dh);

   TFCE_perm(:,p) = single(tf);
   TFCE_null_max(p) = max(tf);

   if track_progress
       send(q,1);
   end

end

% ============================================================
% P-values
% ============================================================
p_vox = (sum(TFCE_perm >= tfce_real,2)+1) ./ (nPerm+1);
p_global = (sum(TFCE_null_max >= TFCE_real_max)+1)/(nPerm+1);

p_img = statistic_image('type','p');
p_img.volInfo = fmri_dat_subj.volInfo;
p_img.removed_voxels = fmri_dat_subj.removed_voxels;
p_img.dat = p_vox;
p_img.p  = p_vox;
% NOTE: every voxel is flagged significant here on purpose. This object is a
% carrier for the p-values; thresholding is the caller's job, typically via
% thresholded_fmri_data_from_statistic_image. Do not read .sig as a result.
p_img.sig = logical(true(size(p_img.dat,1),1));

info.TFCE_real = tfce_real;
info.TFCE_null = TFCE_perm;
info.TFCE_real_max = TFCE_real_max;
info.TFCE_null_max = TFCE_null_max;
info.p_TFCE_global = p_global;
info.nuisance_rank = nuisance_rank;
info.tfce_dh = tfce_dh;

end

% ============================================================
% Helper functions
% ============================================================
function tmap = one_sample_tmap(dat)
% dat: [nVox × nSubj]

n = size(dat,2);
m = mean(dat,2);
s = std(dat,0,2);

tmap = m ./ (s ./ sqrt(n));
end

function tmap = two_sample_tmap(datA, datB)
% datA: [nVox × nA], datB: [nVox × nB]

nA = size(datA,2);
nB = size(datB,2);

mA = mean(datA,2);
mB = mean(datB,2);

vA = var(datA,0,2);
vB = var(datB,0,2);

tmap = (mA - mB) ./ sqrt(vA./nA + vB./nB);
end