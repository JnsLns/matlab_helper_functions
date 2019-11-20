% function mns_out = randnParams(mns,sds)
%
% Draw pseudrandom value from normal distribution with defined mean and
% standard deviation.
%
% Takes column vectors (!) of means and associated standard deviations as
% input. Generates a column vector of same size with values drawn from a
% normal distribution defined by the parameters from the resepective rows
% of the input vectors.
%
% mns = column vector with means
% sds = column vector with stds
% mns_out = column vector with values drawn from normal distr. with input
% parameters
%
% Note:
%
% The general theory of random variables states that
% if x is a random variable with
% mean mu_x and
% variance s^2_x,
%
% then the random variable y,
% defined by y = ax + b,
% where a and b are constants, has
% mean mu_y = a * mu_x + b and
% variance s^2_y = a^2 * s^2_x.


function mns_out = randnParams(mns,sds)

mns_out = sds.*randn(size(sds,1),1) + mns;

end


