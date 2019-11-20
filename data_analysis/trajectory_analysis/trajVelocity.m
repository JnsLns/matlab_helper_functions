function vCell = trajVelocity(trajData, coordCols, tCol, treatmentOfTopRow, customTopRow)
%
% Takes cell array of trajectory data as input, each cell holding one
% trajectory in the form of rows representing timesteps and columns
% representing different position coordinate axes plus time points that
% correspond to the position data.
% For each of these trajectories, this function computes the
% velocity at each time point as the difference between the data
% in one row and that of the preceding one (for the columns provided in
% coordCols and tCol), getting the distance moved and time elapsed between
% these measurements, and then dividing the distance moved by the time
% elapsed between the data points.
%
% The returned cell array contains matrices identical to the ones in the
% input cell array, but with the columns denoted by coordCols replaced by
% the velocity. If treatmentOfTopRow is set to 'drop' (default), the
% matrices will be one row shorter, missing the formerly first row; else the
% first row will be identical to the input except for the columns in
% coordCols, which will be replaced depending on the option chosen (see
% below). If treatmentOfTopRow is 'drop', the first velocity will be paired
% with the very first time point in tCol of the input data (i.e., the
% velocity provided in a row refers to the time interval between the
% time point in the same row and the time point in the next row).
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
% tCol (optional, default 3): scalar value denoting which column in the
% trajectory matrices holds time data.
%
% treatmentOfTopRow (optional, default 'drop'): If set to 'drop' (default),
% the first row will be discarded, output matrices will be one row shorter
% than input. In this case, the velocity provided in a row refers to the
% time interval between the time point in the same row and the time point
% in the next row. If 'zeros', first row in columns given by coordCols will be
% zero; if 'nans', respective elements will be NaN; if 'custom' the last
% argument customTopRow must be provided in the form of an n-element row
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
    tCol = 3;
end

if nargin < 4
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


% Compute moved distance and time intervals btw data points
dxs = trajRowDiffs(trajData, coordCols, 'drop');
dts = trajRowDiffs(trajData, tCol, 'drop');

% divide travelled distances by time (and add desired top row)
dxOverDt = cellfun(@(dx,dt)  [topRow; dx(:,coordCols) ./ repmat(dt(:,tCol),1,numel(coordCols))], dxs,dts,'uniformoutput',0);


% if top row dropped ('drop'), remove it in original data as well.
% (first move t column in that data  down one row, so velocity values will
% be in the same row as the data that represents the starting point of the
% vector for that velocity).
if isempty(topRow)        
    for i = 1:numel(trajData)           
        trajData{i}(:,tCol) = circshift(trajData{i}(:,tCol),1,1);
    end    
    trajData = cellfun(@(tr) tr(2:end,:), trajData, 'uniformoutput', 0);
end

% Write back to original trajectory matrices
for i = 1:numel(dxOverDt)           
    trajData{i}(:,coordCols) = dxOverDt{i};    
end

vCell = trajData;


