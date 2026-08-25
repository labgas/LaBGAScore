function yp = swapWithinSubjectLabels(y, subjIdx)
% swapWithinSubjectLabels  Restricted permutation for a paired design.
%
%   yp = swapWithinSubjectLabels(y, subjIdx)
%
%   Independently and with probability 0.5 swaps the two class labels within
%   each subject. This is the exchangeability-correct null for a within-subject
%   pre/post design: it destroys the condition effect while preserving both the
%   pairing and every subject-level characteristic.
%
%   Because the swap happens WITHIN subject, any subject-constant covariate
%   keeps its exact relationship to the data, so this scheme needs no
%   Freedman-Lane or strata machinery when covariates are supplied.
%
%   Inputs
%     y        [n x 1] binary labels
%     subjIdx  [n x 1] subject index, exactly two rows per subject
%
%   Output
%     yp       [n x 1] permuted labels
%
%   Errors if any subject does not have exactly two rows.
%
%   See also PLSDA_paired_neuroimaging_pipeline, quickGroupedCV.
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
