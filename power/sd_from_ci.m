function SD = sd_from_ci(CI_low, CI_high, n)
%SD_FROM_CI  Compute standard deviation from a 95% CI around the mean.
%
%   SD = SD_FROM_CI(CI_low, CI_high, n)
%
%   Inputs:
%       CI_low  - lower bound of the 95% confidence interval
%       CI_high - upper bound of the 95% confidence interval
%       n       - sample size
%
%   Output:
%       SD      - estimated standard deviation
%
%   Formula:
%       CI = mean ± t * (SD / sqrt(n))
%       => SD = CI_half_width * sqrt(n) / t
%
%   Note: assumes a two-sided 95% CI and uses the t-distribution
%

    % Half-width of the CI
    CI_half_width = (CI_high - CI_low) / 2;

    % Degrees of freedom
    df = n - 1;

    % Critical t-value for 95% CI (two-tailed)
    t_crit = tinv(0.975, df);

    % Standard deviation
    SD = CI_half_width * sqrt(n) / t_crit;
end