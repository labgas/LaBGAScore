function maxLV = capLV(maxLVopt, Xtrain)
% capLV  Cap the number of PLS latent variables to what the data support.
%
%   maxLV = capLV(maxLVopt, Xtrain)
%
%   Returns the largest latent-variable count that is valid for this training
%   fold: min(requested, rank(Xtrain)-1, nTrain-2), floored at 0.
%
%   Inputs
%     maxLVopt  requested maximum (opts.maxLV)
%     Xtrain    [nTr x p] training-fold features, AFTER preprocessing
%
%   Output
%     maxLV     usable latent-variable count; 0 means the fold should be skipped
%
%   Pass the POST-preprocessing matrix. Residualizing on covariates lowers the
%   rank by up to rank(C), so capping on the raw matrix can request more latent
%   variables than the residualized data actually support.
%
%   The rank call is a full SVD, so this is not free; call it once per fold
%   rather than inside the latent-variable loop.
%
%   See also plsregress, foldPreprocess, PLSR_neuroimaging_pipeline.
nTr = size(Xtrain,1);

% rank-based cap is important even when p<n
rX = rank(Xtrain);

% plsregress needs lv <= min(rank(X), n-1) typically; keep your -1/-2 safeguards
maxLV = min([maxLVopt, rX-1, nTr-2]);

if isnan(maxLV) || isinf(maxLV)
    maxLV = 0;
end

maxLV = floor(maxLV);
end
