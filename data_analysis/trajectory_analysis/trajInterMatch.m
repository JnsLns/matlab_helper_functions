% function trajOut = trajInterMatch(trajData,samplingPeriod,coordCols,tCol,method)
%
% Interpolates n-dimensional trajectories each contained in one element of
% a cell array trajData, returning a similar cell array with interpolated
% data.
% The output matrices consist of columns corresponding to the
% columns specified in coordCols (in the same order), plus a column with 
% the query points as last column.
% Interpolation is done independently for each data column given in
% coordCols (referring to the columns of the trajectory matrices in the
% cell array), using as basis the time values given in column tCol. 
% Values are generated for time points separated by samplingPeriod duration,
% for each time point from trajectory onset (first data point) and offset
% (last data point; including these points). All time values in the output
% trajectories are relative to the onset of the respective trajectory. In
% effect, output trajectories will have different numbers of interpolated
% data points, depending on total trajectory duration, but the values of the
% sampling points will be congruent across trajectories (without the need for
% time normalization).
%
% trajData: n-by-m cell array holding in each element one trajectory in the
% form of a p-by-q matrix: In this matrix, p = data points, and q = point
% coordinates. trajData may also be supplied as a p-by-q matrix of that
% form holding only one trajectory.
%
% samplingPeriod: Scalar specifying the temporal interval by which generated
% data points will be separated. Sample points will be tMin:samplingPeriod:tMax,
% where tMin is 0 and tMax is the trajectory duration.
% Note that if mod(tMax,samplingPeriod) ~= 0, the excess data will not be
% represented in the output trajectory.
%
% coordCols (optional, default [1 2 3]): n-element row vector giving the
% column numbers of those columns in each trajectory matrix that hold
% the coordinates that should be interpolated.
%
% tCol (optional, default numel(coordCols)+1): integer giving the column
% number of the column in the trajectory matrices that holds time values.
%
% method (optional, default 'spline'): Interpolation method string. See doc
% interp1.


function trajOut = trajInterMatch(trajData,samplingPeriod,coordCols,tCol,method)

if isa(trajData,'double')            
    trajData = {trajData};            
end

if any(cellfun(@isempty, trajData))
    error('Cell array must not have empty cells')
end

if nargin < 3
    coordCols = [1 2 3];
end

if nargin < 4
    tCol = numel(coordCols)+1;
end

if nargin < 5
    method = 'spline';
end

% ensure that trajData has only one column (avoids problems with cellfun
% outputs, permute, and in case of interaction with other functions)
trajData = reshape(trajData,numel(trajData),1);

% Get min and max t values for each trajectory
tMinMax = cell2mat(cellfun(@(tr)[min(tr(:,tCol)),max(tr(:,tCol))],trajData,'uniformoutput',0));
durations = diff(tMinMax,1,2);

% Generate query points (relative to tZero) for each trajectory
qps = arrayfun(@(x) 0:samplingPeriod:x, durations,'uniformoutput',0);

% Do interpolation (interp is done relative to tZero)
% (note: reshape is required to cover both cases with one and
% multiple columns that are interpolated since interp1 returns a row vector
% in one, and a matrix in the other case)
if ~isempty(tMinMax)
    trajOut = ...
        cellfun(@(tr,qp,tZero) [reshape(interp1(tr(:,tCol)-tZero,tr(:,coordCols),qp,method),[],numel(coordCols)),qp'], ...
        trajData,qps,num2cell(tMinMax(:,1)),'uniformoutput',0);
else
    trajOut = cell([],1); 
end

end

