function results = PLSDA_paired_neuroimaging_pipeline(X, Y, subjectID, opts)
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
%   opts.seed     default 1
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
%     .cvObserved, .cvPredicted, .cvSubjectID, .cvRepeatID
%
% NOTES
%   - This pipeline assumes exactly 2 observations per subject and one
%     observation in each class. It will error otherwise.
%   - This is the correct way to do "raw pre/post row" discrimination
%     without subject leakage across folds.
%
% See also: plsregress, perfcurve

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

if numel(subjectID) ~= n
    error('subjectID must have one entry per row of X.');
end
subjectID = subjectID(:);

[uniqueSubj,~,subjIdx] = unique(subjectID,'stable');
nSubj = numel(uniqueSubj);

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

        Xtrain = X(trainIdx,:);
        Xtest  = X(testIdx,:);

        [Xtrain, Xtest] = applyScaling(Xtrain, Xtest, opts.scale);

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
                [Xtr, Xva] = applyScaling(Xtr, Xva, opts.scale);

                [~,~,~,~,beta] = plsregress(Xtr, ytr, lv);

                yhat = [ones(sum(va),1) Xva] * beta;
                [~,~,~,foldAUC(f)] = perfcurve(yva, yhat, 1);
            end

            innerAUC(lv) = nanmean(foldAUC);
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
results.AUC     = nanmean(AUC(:));

results.allACC  = ACC;
results.ACC     = nanmean(ACC(:));

results.allSENS = SENS;
results.SENS    = nanmean(SENS(:));

results.allSPEC = SPEC;
results.SPEC    = nanmean(SPEC(:));

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

[Xz,~] = applyScaling(X, X, opts.scale);

finalLV = round(nanmedian(selectedLV(:)));
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

permAUC = nan(opts.nPerm,1);

for i = 1:opts.nPerm
    yp = swapWithinSubjectLabels(yNum, subjIdx);
    permAUC(i) = quickGroupedCV(X, yp, subjIdx, outerK, opts);
end

results.allpermAUC = permAUC;
results.permAUC = nanmean(permAUC);
results.permutation_p = mean(permAUC >= results.AUC, 'omitnan');

%% -------------------------------------------------
% 7. Optional simple permutation histogram
%% -------------------------------------------------

figure('Color','w');
histogram(permAUC(~isnan(permAUC)));
hold on;
xline(results.AUC,'LineWidth',1.5);
xlabel('AUC');
ylabel('Frequency');
title('Paired permutation AUC distribution');
grid on;

end

%% =================================================
% Helpers
%% =================================================

function foldID = makeGroupedFolds(nGroups, K)
% Randomly assign groups to K folds as evenly as possible
perm = randperm(nGroups);
foldID = nan(nGroups,1);
for i = 1:nGroups
    foldID(perm(i)) = mod(i-1, K) + 1;
end
end

function maxLV = capLV(maxLVopt, Xtrain)
nTr = size(Xtrain,1);
rX = rank(Xtrain);
maxLV = min([maxLVopt, rX-1, nTr-2]);
if isnan(maxLV) || isinf(maxLV)
    maxLV = 0;
end
maxLV = floor(maxLV);
end

function [XtrS, XteS] = applyScaling(Xtr, Xte, mode)
switch lower(mode)
    case 'none'
        XtrS = Xtr;
        XteS = Xte;
    case 'center'
        mu = mean(Xtr,1);
        XtrS = Xtr - mu;
        XteS = Xte - mu;
    otherwise
        mu = mean(Xtr,1);
        sd = std(Xtr,0,1);
        sd(sd == 0) = 1;
        XtrS = (Xtr - mu) ./ sd;
        XteS = (Xte - mu) ./ sd;
end
end

function yp = swapWithinSubjectLabels(y, subjIdx)
% Under the paired null, randomly swap the two labels within each subject
yp = y;
nSubj = max(subjIdx);

for s = 1:nSubj
    idx = find(subjIdx == s);
    if numel(idx) ~= 2
        error('swapWithinSubjectLabels assumes exactly 2 rows per subject.');
    end

    if rand < 0.5
        yp(idx) = yp(flipud(idx));
    end
end
end

function AUC = quickGroupedCV(X, Y, subjIdx, K, opts)
% Lightweight grouped CV used for paired permutation testing

nSubj = max(subjIdx);
K = min(K, nSubj);
foldID = makeGroupedFolds(nSubj, K);

auc = nan(K,1);

for k = 1:K
    teSubj = foldID == k;
    trSubj = ~teSubj;

    tr = trSubj(subjIdx);
    te = teSubj(subjIdx);

    ytr = Y(tr);
    yte = Y(te);

    if numel(unique(ytr)) < 2 || numel(unique(yte)) < 2
        continue
    end

    Xtr = X(tr,:);
    Xte = X(te,:);
    [Xtr, Xte] = applyScaling(Xtr, Xte, opts.scale);

    lv = capLV(opts.maxLV, Xtr);
    if lv < 1
        continue
    end

    [~,~,~,~,beta] = plsregress(Xtr, ytr, lv);
    score = [ones(sum(te),1) Xte] * beta;
    [~,~,~,auc(k)] = perfcurve(yte, score, 1);
end

AUC = nanmean(auc);
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

