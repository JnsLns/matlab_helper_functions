function [polMat, theta_qps, rho_qps] = cart2pol2dMat(crtn,xVals,yVals,theta_qNum,rho_qNum,method)
%
% Convert 2-dimensional data over cartesian space into representation over
% polar space. The resulting 2d matrix polMat is defined over the interval
% [-pi,pi] along columns 1:end, and over [0,maxRho] along rows 1:end, where
% maxRho is the greatest radius among the coordinates of the supplied
% sample points.
% Optional outputs theta_qps and rho_qps are polar coordinates of the 
% query points used, along the angular and radial dimension, respectively.
% Angle is measured counter-clockwise starting at the positive x-axis
% (i.e., to the right).
%
% crtn: 2d matrix, where columns span x-axis and rows span y-axis. 
%
% xVals and yVals (optional, [] for default): row vectors with a number of
% elements equal to the column and row number of crtn, respectively, 
% holding x and y coordinates of the data in crtn. The values in the two
% vectors may change non-monotonously.
% xVals and yVals also set the location of the origin relative to the
% cartesian data points (which is also the origin w.r.t which the angle and
% radius are determined for the conversion to polar). The origin may lie
% in between the supplied data points.
% If only xVals, only yVals, or neither is supplied, the sampling coordi-
% nates for that dimension are assumed to change linearly and the origin
% is placed in the center of the data along the respective dimension (note
% that if the size of crtn along that dimension is an even number, this
% places the origin between the two middle rows or columns).
%
% theta_qNum and rho_qNum: numnber of evenly spaced query points over angle
% and radius, respectively; determines size of the output matrix along the
% angular dimension (columns) and the radial dimension (rows).  
%
% method (optional, default 'linear'): interpolation method. See
% doc scatteredInterpolant for possible methods.

% size of input matrix
[ySz,xSz] = size(crtn);

if isempty(xVals)            
    xVals = linspace(-xSz/2,xSz/2,xSz);
end

if isempty(yVals)
    yVals = linspace(-ySz/2,ySz/2,ySz);
end

if nargin < 6 
    method = 'linear';
end          

% make matrix same size as crtn, holding x/y values for each corresponding 
% data point in crtn.
[xGrid,yGrid] = meshgrid(xVals,yVals);

% convert these cartesian coordinates to polar coordinates
[theta,rho] = cart2pol(xGrid,yGrid);

% make query points over angle and radius
theta_qps = linspace(-pi,pi,theta_qNum);
rho_qps = linspace(0,max(max(abs(rho))),rho_qNum);

% make interpolant function that can be used to compute query points
% (note: interp2 cannot be used since it requires monotonic sample
% coordinates)
F = scatteredInterpolant(theta(:),rho(:),crtn(:),method);

% Make grid of query points
[theta_qpGrid,rho_qpGrid] = meshgrid(theta_qps,rho_qps);

% Evaluate interpolant for that grid
polMat = F(theta_qpGrid,rho_qpGrid);