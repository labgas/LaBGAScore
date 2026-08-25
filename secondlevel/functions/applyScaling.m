function [XtrS, XteS, sc] = applyScaling(Xtr, Xte, mode)
% applyScaling  Leakage-free feature scaling fit on the training fold.
%
%   [XtrS, XteS]     = applyScaling(Xtr, Xte, mode)
%   [XtrS, XteS, sc] = applyScaling(Xtr, Xte, mode)
%
%   Computes the scaling constants from Xtr alone and applies them to both Xtr
%   and Xte. Called inside every cross-validation fold of the ML pipeline
%   family; scaling the full feature matrix before the CV loop would leak test
%   -fold information and silently invalidate every performance estimate.
%
%   Inputs
%     Xtr   [nTr x p]  training features
%     Xte   [nTe x p]  test features (may be [])
%     mode  'zscore' (default) | 'center' | 'none'
%
%   Outputs
%     XtrS, XteS  scaled features
%     sc          struct with .mode, .mu and .sd (the training-fold constants)
%
%   Zero-variance columns are given sd = 1 so they pass through as zeros rather
%   than producing Inf or NaN.
%
%   See also foldPreprocess, residualizeFold, capLV.

switch lower(mode)
    case 'none'
        XtrS = Xtr;
        XteS = Xte;
        sc = struct('mode','none','mu',[],'sd',[]);

    case 'center'
        mu = mean(Xtr,1);
        XtrS = Xtr - mu;
        if isempty(Xte), XteS = Xte; else, XteS = Xte - mu; end
        sc = struct('mode','center','mu',mu,'sd',[]);

    otherwise % 'zscore'
        mu = mean(Xtr,1);
        sd = std(Xtr,0,1);
        sd(sd==0) = 1;
        XtrS = (Xtr - mu) ./ sd;
        if isempty(Xte), XteS = Xte; else, XteS = (Xte - mu) ./ sd; end
        sc = struct('mode','zscore','mu',mu,'sd',sd);
end

end
