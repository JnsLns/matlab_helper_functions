function [ci_borders, moe, df] = confidence_intervals(x, c_level)
% function [ci_borders, moe] = confidence_intervals(x, c_level)
%
% Compute confidence interval(s) at given confidence level.
%
% x is an n-by-m array with n>1. Confidence intervals at the confidence level
% supplied by c_level (e.g., 0.95) are computed seperately for each column
% of x. ci_borders is a 2-row matrix, where the rows hold the lower and
% upper borders of the confidence intervals, respectively, while the columns
% correspond to those of x. moe is a row vector with as many columns as x,
% holding the corresponding margin of error values.
%
% Note: NaN values in x are completely disregarded in the computation of
% confidence intervals.

if size(x,1) < 2
    error('Input data must have at least two rows.')
end

n = size(x,1) - sum(isnan(x)); % Determine number of cases for each column, not counting NaNs

if any(n<2)
    warning(['At least one column of the input data includes less than two non-NaN cases! ', ...
        'Output will be set to NaN for these.'])
end

se = nanstd(x) ./ sqrt(n);          % standard error for each column
df = n-1;                           % degrees of freedom
Z = abs(tinv((1-c_level)/2, df));   % t-value corresponding to confidence level
moe = Z * se;                       % margin of error (half width of confidence interval)
m = nanmean(x);
ci_borders = [m - moe; m + moe];
moe(n<2) = nan;
ci_borders(:,n<2) = nan;

end