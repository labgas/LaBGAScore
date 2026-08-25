function foldID = makeGroupedFolds(nGroups, K)
% makeGroupedFolds  Assign groups (not rows) to K cross-validation folds.
%
%   foldID = makeGroupedFolds(nGroups, K)
%
%   Randomly distributes nGroups group labels across K folds as evenly as
%   possible. In the paired PLS-DA pipeline the groups are SUBJECTS, so both of
%   a subject's observations always land on the same side of a split and the
%   within-subject pairing cannot leak across the train/test boundary.
%
%   Inputs
%     nGroups  number of groups (subjects)
%     K        number of folds
%
%   Output
%     foldID   [nGroups x 1] fold index per group
%
%   Expand to rows with foldID(subjIdx), where subjIdx maps each row to its
%   group.
%
%   See also PLSDA_paired_neuroimaging_pipeline, quickGroupedCV.
perm = randperm(nGroups);
foldID = nan(nGroups,1);
for i = 1:nGroups
    foldID(perm(i)) = mod(i-1, K) + 1;
end
end
