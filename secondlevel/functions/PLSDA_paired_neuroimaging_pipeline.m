function results = PLSDA_paired_neuroimaging_pipeline(X, Y, subjectID, opts)
% PLSDA_paired_neuroimaging_pipeline  Within-subject PLS-DA for paired designs.
%
%   results = PLSDA_paired_neuroimaging_pipeline(X, Y, subjectID, opts)
%
%   Discriminates two conditions measured in the same subjects (pre/post, on/off)
%   from raw observation rows, with cross-validation folds formed over SUBJECTS so
%   the pairing cannot leak across the train/test boundary. Returns the same field
%   names as PLSDA_neuroimaging_pipeline, so plot_PLSDA_diagnostics_neuroimaging
%   consumes its output unchanged.
%
%   See also README_PLSDA_paired_neuroimaging_pipeline.md for the full guide.
%
% INPUTS
%   X          [n x p] feature matrix
%   Y          [n x 1] binary labels (numeric/logical/categorical/string/cellstr)
%              Example: {'pre','post'} or [0 1]
%   subjectID  [n x 1] subject identifier; each subject must appear exactly twice
%              and contribute one observation to each class
%   opts       struct of options (optional)
%
% OPTIONS
%   opts.outerK   default 5
%   opts.innerK   default 4
%   opts.nRepeats default 50
%   opts.maxLV    default 4
%   opts.nPerm    default 1000
%   opts.scale    default 'zscore'   % 'zscore'|'center'|'none'
%   opts.seed     default 1          % results are reproducible from this seed
%
%   opts.covariates      default []  [n x nCov] numeric nuisance matrix, one ROW
%                        per observation (not per subject). Do NOT include a
%                        column of ones. Nuisance coefficients are estimated on
%                        the TRAINING fold only and applied to both folds.
%                        Because trainIdx is a row mask expanded from a subject
%                        mask, the nuisance fit is automatically subject-blocked.
%                        Leave empty to reproduce the previous behaviour exactly.
%   opts.covariateNames  default {'cov1',...} cellstr, provenance only.
%
% OUTPUTS
%   results struct with fields:
%     .AUC, .ACC, .SENS, .SPEC
%     .allAUC, .allACC, .allSENS, .allSPEC
%     .selectedLV
%     .betaStore
%     .featureWeights
%     .meanFeatureWeight
%     .meanBeta, .sdBeta, .stabilityZ, .signStability
%     .VIP
%     .finalLV, .betaFinal, .finalXLoadings, .finalYLoadings
%     .varExplainedX, .varExplainedY
%     .allpermAUC, .permAUC, .permutation_p
%     .quickCV_observed   observed AUC from quickGroupedCV, i.e. the SAME
%                         estimator that generates the null. permutation_p is
%                         computed against this rather than against .AUC, so
%                         both sides of the test use one estimator, and uses
%                         the (sum(...) + 1) / (nValid + 1) convention.
%     .cvObserved, .cvPredicted, .cvSubjectID, .cvRepeatID
%     .covariateInfo      .used .names .nCov .rank .residualizeY .order
%                         .permScheme ('within-subject-swap')
%     .seed               the RNG seed actually used
%
% NOTES
%   - This pipeline assumes exactly 2 observations per subject and one
%     observation in each class. It will error otherwise.
%   - This is the correct way to do "raw pre/post row" discrimination
%     without subject leakage across folds.
%   - Permutation uses swapWithinSubjectLabels, which swaps the two labels
%     within each subject. That is already the exchangeability-correct null
%     here and needs no Freedman-Lane or strata machinery when covariates are
%     supplied, because a within-subject swap preserves every subject-level
%     characteristic exactly.
%
% COVARIATES IN A PAIRED DESIGN
%   - Subject-CONSTANT covariates (age, sex, genotype) are meaningful and are
%     NOT wiped out. Scaling uses the pooled training-fold mean rather than
%     within-subject centering, so between-subject nuisance variance really is
%     present in X. Being constant within subject, such covariates are also
%     orthogonal to the within-subject condition contrast and cannot remove the
%     effect of interest.
%   - Covariates that VARY within subject are the dangerous case. If such a
%     covariate tracks the condition, residualizing X on it removes exactly the
%     signal under test and drives AUC toward chance. The pipeline warns when a
%     covariate both varies within subject and correlates |r| > 0.3 with the
%     class label.
%
% See also: plsregress, perfcurve, foldPreprocess, swapWithinSubjectLabels

%% -------------------------------------------------
% 0. Defaults
%% -------------------------------------------------

if nargin < 4
    opts = struct;
end

if ~isfield(opts,'outerK');   opts.outerK = 5; end
if ~isfield(opts,'innerK');   opts.innerK = 4; end
if ~isfield(opts,'nRepeats'); opts.nRepeats = 50; end
if ~isfield(opts,'maxLV');    opts.maxLV = 4; end
if ~isfield(opts,'nPerm');    opts.nPerm = 1000; end
if ~isfield(opts,'scale');    opts.scale = 'zscore'; end
if ~isfield(opts,'seed');     opts.seed = 1; end

warnUnknownOptions(opts, { ...
    'outerK','innerK','nRepeats','maxLV','nPerm','scale','seed', ...
    'covariates','covariateNames'}, ...
    'PLSDA_paired_neuroimaging_pipeline');

rng(opts.seed,'twister');

%% -------------------------------------------------
% 1. Input checks / outcome preparation
%% -------------------------------------------------

if iscell(Y) || isstring(Y) || iscategorical(Y)
    Y = grp2idx(Y);
end
Y = Y(:);

if numel(unique(Y)) ~= 2
    error('Y must contain exactly 2 classes for paired pre-post discrimination.');
end

% Convert to 0/1 using larger label as positive class
yNum = double(Y == max(Y));

[n,p] = size(X);

if ~all(isfinite(X(:)))
    error('PLSDA_paired_neuroimaging_pipeline:nonFiniteX', ...
        'X contains NaN or Inf. Handle missing values upstream (impute or drop).');
end

% Covariates are resolved ONCE and then passed explicitly alongside X.
[Cov, covInfo] = validateCovariates(opts, n, 'PLSDA_paired_neuroimaging_pipeline');

% Covariate caveat specific to the paired design.
% Subject-CONSTANT covariates (age, sex, genotype) are meaningful here and are
% NOT wiped out: scaling uses the pooled training-fold mean rather than
% within-subject centering, so between-subject nuisance variance really is
% present in X. Being constant within subject, they are also orthogonal to the
% within-subject condition contrast, so they cannot remove the effect of
% interest.
% Covariates that VARY within subject are the dangerous case: if such a
% covariate tracks the condition, residualizing X on it removes exactly the
% signal being tested and drives AUC toward chance.

if numel(subjectID) ~= n
    error('subjectID must have one entry per row of X.');
end
subjectID = subjectID(:);

[uniqueSubj,~,subjIdx] = unique(subjectID,'stable');
nSubj = numel(uniqueSubj);

% Warn about covariates that vary within subject and track the class label:
% residualizing X on those removes the effect of interest.
if ~isempty(Cov)
    for jc = 1:size(Cov,2)
        withinRange = accumarray(subjIdx, Cov(:,jc), [], @(v) max(v)-min(v));
        if any(withinRange > 0)
            rho = corr(Cov(:,jc), yNum);
            if abs(rho) > 0.3
                warning('PLSDA_paired_neuroimaging_pipeline:covariateTracksLabel', ...
                    ['Covariate "%s" varies within subject and correlates r = %.2f with the class label. ' ...
                     'Residualizing X on it will remove signal of interest.'], covInfo.names{jc}, rho);
            end
        end
    end
end

if nSubj < 2
    error('Need at least 2 subjects.');
end

% Validate strict paired structure
for s = 1:nSubj
    idx = subjIdx == s;

    if sum(idx) ~= 2
        error('Each subject must appear exactly twice. Subject %s has %d rows.', ...
            toString(uniqueSubj(s)), sum(idx));
    end

    ys = yNum(idx);
    if numel(unique(ys)) ~= 2
        error(['Each subject must contribute one observation to each class. ', ...
               'Subject %s violates this.'], toString(uniqueSubj(s)));
    end
end

% Effective K cannot exceed number of subjects
outerK = min(opts.outerK, nSubj);
innerK_default = min(opts.innerK, max(2, nSubj - ceil(nSubj/outerK)));

%% -------------------------------------------------
% 2. Repeated nested grouped CV
%% -------------------------------------------------

AUC  = nan(opts.nRepeats, outerK);
ACC  = nan(opts.nRepeats, outerK);
SENS = nan(opts.nRepeats, outerK);
SPEC = nan(opts.nRepeats, outerK);

selectedLV    = nan(opts.nRepeats, outerK);
betaStore     = nan(p+1, outerK, opts.nRepeats);
featureWeights = nan(p, opts.nRepeats * outerK);

cvObserved = [];
cvPredicted = [];
cvSubjectID = [];
cvRepeatID = [];

for r = 1:opts.nRepeats

    outerFoldID = makeGroupedFolds(nSubj, outerK);

    for k = 1:outerK

        testSubj  = outerFoldID == k;
        trainSubj = ~testSubj;

        trainIdx = trainSubj(subjIdx);
        testIdx  = testSubj(subjIdx);

        ytrain = yNum(trainIdx);
        ytest  = yNum(testIdx);

        if numel(unique(ytrain)) < 2 || numel(unique(ytest)) < 2
            continue
        end

        Xtrain_raw = X(trainIdx,:);
        Xtest_raw  = X(testIdx,:);

        Ctrain = []; Ctest = [];
        if ~isempty(Cov)
            Ctrain = Cov(trainIdx,:);
            Ctest  = Cov(testIdx,:);
        end

        % trainIdx is a row logical expanded from a subject mask, so slicing Cov
        % with it keeps the nuisance fit subject-blocked for free.
        [Xtrain, Xtest] = foldPreprocess(Xtrain_raw, Xtest_raw, Ctrain, Ctest, opts.scale);

        maxLV = capLV(opts.maxLV, Xtrain);
        if maxLV < 1
            continue
        end

        % ----- inner grouped CV on training subjects only -----
        trainSubjList = find(trainSubj);
        nTrainSubj = numel(trainSubjList);
        innerK = min(innerK_default, nTrainSubj);
        innerK = max(2, innerK);

        innerFoldID_local = makeGroupedFolds(nTrainSubj, innerK);
        innerFoldID_global = nan(nSubj,1);
        innerFoldID_global(trainSubjList) = innerFoldID_local;

        innerAUC = nan(maxLV,1);

        for lv = 1:maxLV
            foldAUC = nan(innerK,1);

            for f = 1:innerK
                vaSubj = innerFoldID_global == f;
                trSubj = trainSubj & ~vaSubj;

                tr = trSubj(subjIdx);
                va = vaSubj(subjIdx);

                ytr = yNum(tr);
                yva = yNum(va);

                if numel(unique(ytr)) < 2 || numel(unique(yva)) < 2
                    continue
                end

                Xtr = X(tr,:);
                Xva = X(va,:);
                C2 = []; Cv = [];
                if ~isempty(Cov), C2 = Cov(tr,:); Cv = Cov(va,:); end
                [Xtr, Xva] = foldPreprocess(Xtr, Xva, C2, Cv, opts.scale);

                [~,~,~,~,beta] = plsregress(Xtr, ytr, lv);

                yhat = [ones(sum(va),1) Xva] * beta;
                [~,~,~,foldAUC(f)] = perfcurve(yva, yhat, 1);
            end

            innerAUC(lv) = mean(foldAUC,'omitnan');
        end

        [~,bestLV] = max(innerAUC);
        selectedLV(r,k) = bestLV;

        % ----- final fit on outer training set -----
        [~,~,~,~,beta] = plsregress(Xtrain, ytrain, bestLV);

        betaStore(:,k,r) = beta;
        featureWeights(:,(r-1)*outerK + k) = beta(2:end);

        scores = [ones(sum(testIdx),1) Xtest] * beta;

        [~,~,~,AUC(r,k)] = perfcurve(ytest, scores, 1);
        pred = scores > 0.5;

        ACC(r,k) = mean(pred == ytest);

        tp = sum(pred == 1 & ytest == 1);
        tn = sum(pred == 0 & ytest == 0);
        fp = sum(pred == 1 & ytest == 0);
        fn = sum(pred == 0 & ytest == 1);

        if (tp + fn) > 0
            SENS(r,k) = tp / (tp + fn);
        end
        if (tn + fp) > 0
            SPEC(r,k) = tn / (tn + fp);
        end

        cvObserved = [cvObserved; ytest]; %#ok<AGROW>
        cvPredicted = [cvPredicted; scores]; %#ok<AGROW>
        cvSubjectID = [cvSubjectID; subjectID(testIdx)]; %#ok<AGROW>
        cvRepeatID = [cvRepeatID; repmat(r, sum(testIdx), 1)]; %#ok<AGROW>
    end
end

results.allAUC  = AUC;
results.AUC     = mean(AUC(:),'omitnan');

results.allACC  = ACC;
results.ACC     = mean(ACC(:),'omitnan');

results.allSENS = SENS;
results.SENS    = mean(SENS(:),'omitnan');

results.allSPEC = SPEC;
results.SPEC    = mean(SPEC(:),'omitnan');

results.selectedLV = selectedLV;
results.betaStore = betaStore;
results.featureWeights = featureWeights;
results.meanFeatureWeight = mean(featureWeights, 2, 'omitnan');

results.cvObserved = cvObserved;
results.cvPredicted = cvPredicted;
results.cvSubjectID = cvSubjectID;
results.cvRepeatID = cvRepeatID;

fprintf('Grouped paired nested CV AUC = %.3f\n', results.AUC);

%% -------------------------------------------------
% 3. Final model (interpretation only)
%% -------------------------------------------------

% Interpretation-only fit on all rows; no held-out set, so a full-sample
% nuisance fit is correct by construction.
[Xz,~] = foldPreprocess(X, [], Cov, [], opts.scale);

finalLV = round(median(selectedLV(:),'omitnan'));
finalLV = max(1, min(finalLV, capLV(opts.maxLV, Xz)));

[XL_final,YL_final,XS,~,betaFinal,PCTVAR,~,stats] = plsregress(Xz, yNum, finalLV);

results.finalLV = finalLV;
results.betaFinal = betaFinal;
results.varExplainedX = PCTVAR(1,:);
results.varExplainedY = PCTVAR(2,:);
results.finalXLoadings = XL_final;
results.finalYLoadings = YL_final;

%% -------------------------------------------------
% 4. VIP
%% -------------------------------------------------

W = stats.W;
T = XS;
Q = YL_final;

SSY = sum(T.^2,1) .* (Q'.^2);

VIP = zeros(p,1);
for j = 1:p
    w = (W(j,:).^2) ./ sum(W.^2,1);
    VIP(j) = sqrt(p * sum(SSY .* w) / sum(SSY));
end

results.VIP = VIP;

%% -------------------------------------------------
% 5. Stability metrics
%% -------------------------------------------------

betaMat = reshape(betaStore(2:end,:,:), p, []);
meanBeta = mean(betaMat, 2, 'omitnan');
sdBeta   = std(betaMat, 0, 2, 'omitnan');
sdBeta(sdBeta == 0) = NaN;

results.meanBeta = meanBeta;
results.sdBeta = sdBeta;
results.stabilityZ = meanBeta ./ sdBeta;

meanSign = sign(results.meanFeatureWeight);
meanSign(meanSign == 0) = NaN;

betaSign = sign(betaMat);
betaSign(betaSign == 0) = NaN;

results.signStability = mean(meanSign == betaSign, 2, 'omitnan');

%% -------------------------------------------------
% 6. Paired permutation test (within-subject label swap)
%% -------------------------------------------------

% Matched observed statistic: the null is generated by quickGroupedCV, so the
% value it is compared against must come from quickGroupedCV too.
obsMatched = quickGroupedCV(X, yNum, subjIdx, Cov, outerK, opts);
results.quickCV_observed = obsMatched;

permAUC = nan(opts.nPerm,1);

for i = 1:opts.nPerm
    yp = swapWithinSubjectLabels(yNum, subjIdx);
    permAUC(i) = quickGroupedCV(X, yp, subjIdx, Cov, outerK, opts);
end

results.allpermAUC = permAUC;
results.permAUC = mean(permAUC,'omitnan');
results.permutation_p = (sum(permAUC >= obsMatched) + 1) / (sum(~isnan(permAUC)) + 1);

%% -------------------------------------------------
% 7. Optional simple permutation histogram
%% -------------------------------------------------

figure('Color','w');
histogram(permAUC(~isnan(permAUC)));
hold on;
xline(obsMatched,'LineWidth',1.5);
xlabel('AUC');
ylabel('Frequency');
title('Paired permutation AUC distribution');
grid on;

%% =================================================
% Provenance
%% =================================================

covInfo.residualizeY = false;   % never applicable: Y is binary here
covInfo.order        = 'residualize-then-scale';
covInfo.permScheme   = 'within-subject-swap';
results.covariateInfo = covInfo;
results.seed          = opts.seed;

end

function s = toString(x)
if isstring(x) || ischar(x)
    s = char(x);
elseif isnumeric(x) || islogical(x)
    s = num2str(x);
else
    s = '<subject>';
end

end
