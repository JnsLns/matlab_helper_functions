function diffCell = trajRowDiffs(trajData, coordCols, treatmentOfTopRow, customTopRow)
%
% Takes cell array of trajectory data as input, each cell holding one
% trajectory in the form of rows representing timesteps and columns
% representing different coordinate axes.
% For each of these trajectories, this function computes the difference be-
% tween the data in one row and that of the preceding one (only for the
% columns provided in coordCols). 
%
% The returned cell array contains matrices identical to the ones in the
% input cell array, but with the columns denoted by coordCols replaced by
% the row differences. If treatmentOfTopRow is set to 'drop' (default), the
% matrices will be one row shorter, missing the formerly first row; else the
% first row will be identical to the input except for the columns in
% coordCols, which will be replaced depending on the option chosen (see
% below).
%
% trajData: n-by-m cell array holding in each element one 2D trajectory in
% the form of an n-by-m matrix: In this matrix, n = data points, and
% m = point coordinates or other data types (column numbers of actual
% x/y coordinates are provided in coordCols). trajData may also be supplied
% as an n-by-m matrix of that form holding only one trajectory (the function
% still returns a cell array).
%
% coordCols (optional, default [1 2]): n-element row vector giving the
% column numbers of the columns for which differences should be calculated.
%
% treatmentOfTopRow (optional, default 'drop'): If set to 'drop' (default),
% the first row will be discarded, output matrices will be one row shorter
% than input. If 'zeros', first row in columns given by coordCols will be
% zero; if 'nans', respective elements will be NaN; if 'custom' the last
% argument customTopRow must be provided in the form of a n-element row
% vector where n == numel(coordCols), the first row will be replaced by that
% vector.

if isa(trajData,'double')
    trajData = {trajData};
end

if any(cellfun(@isempty, trajData))
    error('Cell array must not have empty cells')
end

if nargin < 2 || all(coordCols == 0)
    coordCols = [1 2];
end

if nargin < 3
    treatmentOfTopRow = 'drop';
end

% Fill-in for first row 
switch treatmentOfTopRow
    case 'drop'        
        topRow = [];
    case 'zeros'
        topRow = zeros(1,numel(coordCols));
    case 'nans'
        topRow = nan(1,numel(coordCols));
    case 'custom'
        topRow = customTopRow;
end
    
% get relevant columns from data
tmpDat = cellfun(@(tr)  tr(:,coordCols) ,trajData,'uniformoutput',0);

% from each row subtract preceding one and add top row (or not)
tmpDat = cellfun(@(tr)  [topRow; tr(2:end,:) - tr(1:end-1,:)],tmpDat,'uniformoutput',0);

% if top row dropped ('drop'), remove it in original data as well.
if isempty(topRow)
    trajData = cellfun(@(tr) tr(2:end,:), trajData, 'uniformoutput', 0);
end

% Write back to original trajectory matrices
for i = 1:numel(tmpDat)           
    trajData{i}(:,coordCols) = tmpDat{i};    
end

diffCell = trajData;



