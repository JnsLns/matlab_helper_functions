% function trajOut = trajInter(trajData,queryPoints,coordCols,tCol,method)
%
% Interpolates n-dimensional trajectories each contained in one element of
% a cell array trajData, returning a similar cell array with interpolated
% data. The output matrices consist of columns corresponding to the
% columns specified in coordCols (in the same order), plus a column with 
% the query points as last column. 
% Interpolation is done independently for each data column that is pro-
% vided in coordCols. The time values given in column tCol are used as
% sample points. (In other words, data along each dimension is treated as
% a 1-dimensional function of time and interpolated independently from the
% other dimensions using interp1.)
% The function generates values for a number of queryPoints data points.
% These are are equally distributed between (and including) the first and
% last time value in tCol, using the interpolation method denoted by 
% method. The matrices in the ouptut cell each have queryPoints rows.
%
% trajData: n-by-m cell array holding in each element one trajectory in the
% form of an p-by-q matrix: In this matrix, p = data points, and q = point
% coordinates. trajData may also be supplied as an p-by-q matrix of that
% form holding only one trajectory.
%
% queryPoints: *Either* scalar giving the number of desired query points
% (i.e., the number of data points each output trajectory will have). The
% actual query points will be evenly distributed between the first and last
% value in the column specified by tCol. *Or* n-by-m matrix where n is the
% number of trajectories in trajData and m is the number of query points.
% Each row in queryPoints is then used for the corresponding trajectory.
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
%
%

function trajOut = trajInter(trajData,queryPoints,coordCols,tCol,method)

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

% generate query points
if numel(queryPoints) == 1
    qps =...
        cellfun(@(tr)  linspace(tr(1,tCol),tr(end,tCol),queryPoints),...
        trajData,'uniformoutput',0);
elseif numel(queryPoints) > 1
    qps = num2cell(queryPoints,2);
end

% interpolate (note: reshape is required to cover both cases with one and
% multiple columns that are interpolated since interp1 returns a row vector
% in one, and a matrix in the other case)
trajOut = ...
    cellfun(@(tr,qp) [reshape(interp1(tr(:,tCol),tr(:,coordCols),qp,method),[],numel(coordCols)),qp'], ...
    trajData,qps,'uniformoutput',0);

end

