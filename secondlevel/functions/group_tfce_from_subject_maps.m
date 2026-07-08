function [tfce_dat, p_img, info] = group_tfce_from_subject_maps( ...
   fmri_dat_subj, design, group, covariates, nPerm, varargin)
% GROUP_TFCE_FROM_SUBJECT_MAPS
%
% Group-level TFCE permutation inference starting from subject-level
% contrast images, with optional Freedman–Lane covariate control.
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
%
% nPerm : scalar integer
%  - number of permutations
%
% OPTIONAL NAME–VALUE PAIRS
%  'H'       : TFCE height exponent (default = 2)
%  'E'       : TFCE extent exponent (default = 0.5)
%  'conn'    : connectivity (default = 26)
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

voxsize = mean(abs(diag(fmri_dat_subj.volInfo.mat(1:3,1:3))));
assert(isfinite(voxsize) && isscalar(voxsize));

% ============================================================
% Step 1: Residualise subject maps (Freedman–Lane)
% ============================================================
Y = fmri_dat_subj.dat; % [nVox × nSubj]

if use_covariates
   Xn = [ones(nSubj,1) covariates];
   beta_nuis = Xn \ Y';
   Y_hat = (Xn * beta_nuis)';    % fitted nuisance effect
   R = Y - Y_hat;                % residuals
else
   R = Y;
   Y_hat = zeros(size(Y));
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

tfce_real = tfce_one_fmri_dat( ...
   t_real,...
   fmri_dat_subj.removed_voxels,...
   fmri_dat_subj.volInfo,...
   voxsize,H,E,conn,sidedness,tail);

TFCE_real_max = max(tfce_real);

tfce_dat = fmri_dat_subj;
tfce_dat.dat = single(tfce_real);
tfce_dat.dat_descrip = sprintf('TFCE (%s, FL)',design);

% ============================================================
% PERMUTATIONS (Freedman–Lane)
% ============================================================
TFCE_perm = zeros(nVox,nPerm,'single');
TFCE_null_max = zeros(nPerm,1,'single');

if use_parallel && nPerm>1
   pool = gcp('nocreate');
   if isempty(pool), parpool; end
   tracker = ProgressTracker(nPerm);
   q = parallel.pool.DataQueue;
   afterEach(q,@(~) tracker.update());
end

parfor p = 1:nPerm

   % --- permute residuals ---
   switch design
       case 'onesample'
           signs = (rand(nSubj,1)>0.5)*2-1;
           R_perm = R .* signs';

       case 'twosample'
           perm_idx = randperm(nSubj);
           R_perm = R(:,perm_idx);
           gp = group(perm_idx);
           gA = gp == grp_vals(1);
           gB = gp == grp_vals(2);
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
       voxsize,H,E,conn,sidedness,tail);

   TFCE_perm(:,p) = single(tf);
   TFCE_null_max(p) = max(tf);
   
   send(q,1);
   
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
p_img.sig = logical(true(size(p_img.dat,1),1));

info.TFCE_real = tfce_real;
info.TFCE_null = TFCE_null;
info.TFCE_real_max = TFCE_real_max;
info.TFCE_null_max = TFCE_null_max;
info.p_TFCE_global = p_global;

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