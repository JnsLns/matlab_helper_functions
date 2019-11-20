% function [MDs, rowNos, xy_md, xy_dpfs, tValues] = ...
%     trajMaxDev(trajData, xyCols, tCol, nQueryPoints, samplingPeriod, interpMethod)
%
% Computes maximum deviations (MD) for each 2d trajectory in a cell array
% trajData. MD is defined as the highest euclidean distance of any data
% point from a straight line passing through the start and end point of a
% trajectory. The first data point is the one given in the first row of
% each trajectory matrix, the last point is the one in the last row (true
% story). If there are multiple highest data points in a trajectory, only
% the first one is considered for all outputs. All outputs are column
% vectors. Output MDs contains the maximum deviation for each trajectory
% (negative values are deviations to the left of the ideal path from first
% to last point; if in doubt, this should be tested, as it also depends on
% the initial reference frame and axis directionality), tValues contains
% the time points for each value in the other ouputs (in same time frame as
% input data), rowNos provides the row number (i.e. timestep) for each value
% in the other outputs, xy_md contains the corresponding x- and y-Values
% in the same coordinate frame as the input trajectories. xy_dpfs is an
% 2-by-n matrix that holds in each row the coordinates of the dropped
% perpendicular foot from the trajectory point corresponding to max
% deviation onto the straight path from the trajectories first to its last
% data point (useful for plotting this perpendicular).
% By default, all computations are done using non-interpolated versions of
% the input trajectories (the input may already be interpolated though).
% This can be changed by setting nQueryPoints or samplingPeriod (see
% below), to enable interpolation with fixed or variable numbers of query
% points, respectively.
%
% trajData: n-by-m cell array holding in each element one trajectory in the
% form of an p-by-q matrix: In this matrix, p = data points, and q = point
% coordinates. trajData may also be supplied as an p-by-q matrix of that
% form holding only one trajectory.
%
% xyCols (optional, default [1 2]): 2-element row vector giving the
% column numbers of those columns in each trajectory matrix that hold
% the x and y coordinates of the trajectory.
%
% tCol (optional, default 3): integer giving the column number of the
% column in the trajectory matrices that holds time values (required for
% interpolation).
%
% nQueryPoints (optional, default 0): For interplation. Interpolation will
% be done using the function trajInter, so that is results in each
% trajectory having the same number of data points, namely nQueryPoints.
% nQueryPoints can be *either* a scalar giving the number of desired
% query points *or* an n-by-m matrix where n is the number of trajectories
% in trajData and m is the number of query points that will be used for the
% trajectory in the corresponding element. ***If both nQueryPoints and
% samplingPeriod are set to 0, no interpolation will be done before
% calculating outputs!***
%
% samplingPeriod (optional, default 0): For interpolation. If 0, does
% nothing. If anything else, all settings regarding nQueryPoints will be
% overridden, and interpolation will instead be done using trajInterMatch,
% so that the query points for all trajectories will be at the same time
% points relative to trajectory onset, namely 0:samplingPeriod:lastTimePoint.
% This means that trajectories of differing duration will have different
% numbers of query points (and that some data at the end of trajectories
% may be lost in case mod(totalDuration,samplingPeriod)~=0.
%
% interpMethod (optional, default 'spline'): Method used for interpolation.
% See doc interp1.
%


function [MDs, rowNos, xy_md, xy_dpfs, tValues] = ...
    trajMaxDev(trajData, xyCols, tCol, nQueryPoints, samplingPeriod, interpMethod)

if isa(trajData,'double')
    trajData = {trajData};
end

if nargin < 6
    interpMethod = 'spline';
end

if nargin < 5
    samplingPeriod = 0;
end

if nargin < 4
    nQueryPoints = 0;
end

if nargin < 3
    tCol = 3;
end

if nargin < 2
    xyCols = [1 2];
end

if any(cellfun(@isempty, trajData))
    error('Cell array must not have empty cells')
end

% ensure that trajData has only one column (avoids problems with cellfun
% outputs, permute, and in case of interaction with other functions)
trajData = reshape(trajData,numel(trajData),1);

% Interpolate if desired
if samplingPeriod ~= 0
    trajData = trajInterMatch(trajData,samplingPeriod,xyCols,tCol,interpMethod);
elseif nQueryPoints ~= 0
    trajData = trajInter(trajData,nQueryPoints,xyCols,tCol,interpMethod);
end

% Make parallel to y-axis, shift to origin
trajData_originalFrame = trajData;
[trajData, trajStarts, revRefDirs] = trajRot(trajData,2,xyCols,0,1);

% Get max deviations, corresponding t values, and timesteps
[~,maxNdcs] = cellfun(@(x) max(abs(x(:,xyCols(1)))) ,trajData);

if ~isempty(maxNdcs)
    
    % get max deviations and corresponding time values (latter only if desired)
    if nargout > 4 
        MD_t = cell2mat(cellfun(@(tr,ndx)  tr(ndx,[xyCols(1),tCol]) ,trajData,num2cell(maxNdcs),'uniformoutput',0));
        tValues = MD_t(:,2);
    else
        MD_t = cell2mat(cellfun(@(tr,ndx)  tr(ndx,xyCols(1)) ,trajData,num2cell(maxNdcs),'uniformoutput',0));
    end
    MDs = MD_t(:,1);    
    % Get corresponding xy-values in non-rotated/shifted frame
    if nargout >= 4
        xy_md  = cellfun(@(tr,ndx)  tr(ndx,xyCols) ,trajData_originalFrame,num2cell(maxNdcs),'uniformoutput',0);
    end
    
    % Output row number for max dev data point
    rowNos = maxNdcs;
    
    if nargout >= 5
        % Get dropped perpendicular foot of MD trajectory point onto ideal path
        % in rectified frame (trajectories on axis 2 and starting in origin)
        xy_dpfs = [zeros(numel(trajData),1),cell2mat(cellfun(@(tr,ndx) tr(ndx,xyCols(2)),trajData,num2cell(maxNdcs),'uniformoutput',0))];
        % transform back to original coordinate frame
        xy_dpfs = trajRevRot(xy_dpfs,2,xyCols,revRefDirs,trajStarts);
        
    end
    
else
    [MDs, tValues, rowNos, xy_md, xy_dpfs] = deal([]);
end

end