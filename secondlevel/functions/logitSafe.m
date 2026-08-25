function z = logitSafe(p)
% logitSafe  Logit transform guarded against 0 and 1.
%
%   z = logitSafe(p)
%
%   Clamps p away from the open-interval endpoints before taking
%   log(p./(1-p)), so probabilities of exactly 0 or 1 give large finite values
%   rather than -Inf or +Inf. Used when converting Elastic Net predicted
%   probabilities into scores for perfcurve.
%
%   See also quickCV_ENet, perfcurve.
p = min(max(p, 1e-6), 1-1e-6);
z = log(p/(1-p));
