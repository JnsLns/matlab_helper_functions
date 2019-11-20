function [crtnMat, x_qps, y_qps] = pol2cart2dMat(pol,thetaVals,rhoVals,x_qNum,y_qNum,method)
%
% Convert 2-dimensional data over polar coordinates into representation in
% over 2d cartesian space. The data in the resulting 2d matrix crtnMat is
% defined over the interval [-maxRho,maxRho] along both axes (i.e., x and y,
% corresponding to columns and rows of crtnMat, respectively), where maxRho
% is the greatest (absolute) radius in the supplied polar data.
% Optional outputs x_qps and y_qps are cartesian coordinates of the query
% points used, along the angular and radial dimension, respectively.
%
% pol: 2d matrix whose columns span angular dimension (theta) and whose
% rows span radial dimension (rho).
%
% thetaVals and rhoVals (optional, [] for default): row vectors with a
% number of elements equal to the column and row number of pol, respectively, 
% holding angular and radial coordinates of the data in pol. The values in
% the two vectors may change non-monotonously.
% If only thetaVals, only rhoVals, or neither is supplied, the sample
% coordinates for each missing dimension are assumed to change linearly and
% along row or column numbers, with the origin being assumed in the first
% row of input matrix pol and in the center column, or, if column number is
% even, in between the data of the center columns.
%
% theta_qNum and rho_qNum: number of evenly spaced query points over
% angular and radial dimension, respectively; determines size of the output
% matrix along the angular (columns) and radial dimension (rows).  
%
% method (optional, default 'linear'): interpolation method. See
% doc scatteredInterpolant for possible methods.

% size of input matrix
[rhoSz,thetaSz] = size(pol);

if isempty(thetaVals)                
    thetaVals = linspace(-pi,pi,thetaSz);    
end

if isempty(rhoVals)    
    rhoVals = 0:rhoSz-1;
end

if nargin < 6 
    method = 'linear';
end          

% make matrix same size as pol, holding theta/rho values for each corres-
% ponding data point in pol.
[thetaGrid,rhoGrid] = meshgrid(thetaVals,rhoVals);

% convert these polar coordinates to cartesian coordinates
[x,y] = pol2cart(thetaGrid,rhoGrid);

% make query points over x and y
x_qps = linspace(-max(max(abs(rhoGrid))),max(max(abs(rhoGrid))),x_qNum);
y_qps = linspace(-max(max(abs(rhoGrid))),max(max(abs(rhoGrid))),y_qNum);

% make interpolant function that can be used to compute query points
% (note: interp2 cannot be used since it requires monotonic sample
% coordinates)
F = scatteredInterpolant(x(:),y(:),pol(:),method);

% Make grid of query points
[x_qpGrid,y_qpGrid] = meshgrid(x_qps,y_qps);

% Evaluate interpolant for that grid
crtnMat = F(x_qpGrid,y_qpGrid);