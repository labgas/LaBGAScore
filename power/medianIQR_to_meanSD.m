function [mu, sigma] = medianIQR_to_meanSD(medianVal, IQR, n)
% medianIQR_to_meanSD
%
% Approximate mean and standard deviation from median and IQR.
%
% INPUTS
%   medianVal : median
%   IQR       : interquartile range (Q3 - Q1)
%   n         : sample size (optional but recommended)
%
% OUTPUTS
%   mu        : estimated mean
%   sigma     : estimated standard deviation
%
% METHODS
%   - Mean ≈ median (robust assumption for symmetric distributions)
%   - SD approximation depends on n:
%       if n >= 25:   SD ≈ IQR / 1.35
%       if n < 25:    SD ≈ IQR / (2 * norminv((0.75*n - 0.125)/(n + 0.25)))
%
% REFERENCES
%   Wan et al. (2014), BMC Med Res Methodol
%   Luo et al. (2018), Stat Methods Med Res

if nargin < 3 || isempty(n)
    n = NaN;
end

% --------------------------
% Mean approximation
% --------------------------
mu = medianVal;

% --------------------------
% SD approximation
% --------------------------
if isnan(n) || n >= 25
    % Large sample approximation
    sigma = IQR / 1.35;
else
    % Small sample correction (Wan et al.)
    p = (0.75*n - 0.125) / (n + 0.25);
    z = norminv(p);
    sigma = IQR / (2 * z);
end

end