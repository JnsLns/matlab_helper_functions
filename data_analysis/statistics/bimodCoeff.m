function [b] = bimodCoeff(x) 
%m3 = skew 
%m4 = kurt 
%n = data size 

% Sarle's bimodality coefficient, as computed in Freeman & Dale (2013, "Assessing
% bimodality to detect the presence of a dual cognitive process") and
% corrected by Pfister et al.(2013, "Good things peak in pairs: a note on
% the bimodality coeffecient). 
%
% See also
m3 = skewness(x,0); 
m4 = kurtosis(x,0)-3; % note that -3 is to compute *excess* kurtosis (i.e., difference to k. of normal distr.)
n = length(x); 
b=(m3^2+1)  /  (m4 + 3 * ( (n-1)^2  /  ((n-2)*(n-3))  )); 