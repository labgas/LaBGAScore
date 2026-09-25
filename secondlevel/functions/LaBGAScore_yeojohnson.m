function [yt, lambda, info] = LaBGAScore_yeojohnson(y, lambda)
% LaBGAScore_yeojohnson  Yeo-Johnson power transform, lambda by maximum likelihood.
%
% Normalises a skewed variable that contains NEGATIVE values, which Box-Cox
% cannot: Box-Cox needs strictly positive input, so applying it here would mean
% adding an arbitrary shift constant, and the shift changes the fitted lambda.
% Yeo-Johnson handles both signs natively and reduces to Box-Cox on (y+1) for
% non-negative y.
%
% :Usage:
% ::
%     [yt, lambda] = LaBGAScore_yeojohnson(y);          % lambda by MLE
%     yt           = LaBGAScore_yeojohnson(y, lambda);  % apply a known lambda
%
% :Inputs:
%   **y:** numeric vector, may contain negatives and NaN (NaN passed through).
%
% :Optional Inputs:
%   **lambda:** scalar. If omitted or empty, estimated by profile likelihood
%               over [-5, 5].
%
% :Outputs:
%   **yt:**     transformed vector, standardised to mean 0 / sd 1 so the scale
%               is comparable with the untransformed factor scores.
%   **lambda:** the lambda used.
%   **info:**   struct with .skew_before, .skew_after, .n, .loglik.
%
% Transform:
%   y >= 0, lambda ~= 0 : ((y+1)^lambda - 1) / lambda
%   y >= 0, lambda == 0 : log(y+1)
%   y <  0, lambda ~= 2 : -(((-y+1)^(2-lambda)) - 1) / (2-lambda)
%   y <  0, lambda == 2 : -log(-y+1)
%
% ..
%     Copyright (C) 2026 Lukas Van Oudenhove. GPLv3.
% ..

y  = y(:);
ok = ~isnan(y);
if nargin < 2 || isempty(lambda)
    nll    = @(L) -local_loglik(y(ok), L);
    lambda = fminbnd(nll, -5, 5, optimset('TolX', 1e-6));
end

yt      = nan(size(y));
yt(ok)  = local_transform(y(ok), lambda);

info = struct();
info.n           = sum(ok);
info.loglik      = local_loglik(y(ok), lambda);
info.skew_before = local_skew(y(ok));
info.skew_after  = local_skew(yt(ok));

% Standardise so the transformed factor is on a comparable scale to the
% untransformed one. Monotone, so it changes no rank and no p-value.
mu = mean(yt(ok)); sg = std(yt(ok));
if sg > 0, yt(ok) = (yt(ok) - mu) / sg; end

end % main


function z = local_transform(y, L)
z   = zeros(size(y));
pos = y >= 0;
if abs(L) > eps
    z(pos) = ((y(pos) + 1).^L - 1) / L;
else
    z(pos) = log(y(pos) + 1);
end
if abs(L - 2) > eps
    z(~pos) = -(((-y(~pos) + 1).^(2 - L)) - 1) / (2 - L);
else
    z(~pos) = -log(-y(~pos) + 1);
end
end


function ll = local_loglik(y, L)
% Profile log-likelihood: the Jacobian term is (L-1)*sum(sign(y).*log(|y|+1)).
n  = numel(y);
z  = local_transform(y, L);
s2 = var(z, 1);
if s2 <= 0, ll = -Inf; return, end
ll = -n/2 * log(s2) + (L - 1) * sum(sign(y) .* log(abs(y) + 1));
end


function s = local_skew(x)
m = mean(x); sd = std(x);
if sd == 0, s = 0; else, s = mean((x - m).^3) / sd^3; end
end
