% function trajOut = trajMirror(trajData, coordCols)
%
% Mirrors 2-dimensional trajectories across their main direction (straight
% line from first to last data point; first and last data point remain
% unchanged).
%
% trajData: n-by-m cell array holding in each element one 2D trajectory in
% the form of an n-by-m matrix: In this matrix, n = data points, and
% m = point coordinates or other data types (column numbers of actual
% x/y coordinates are provided in coordCols). trajData may also be supplied
% as an n-by-m matrix of that form holding only one trajectory (the function
% still returns a cell array).
%
% coordCols (optional, default [1 2]): 2-element row vector giving the
% column numbers of the columns in each trajectory matrix that hold [x y]
% coordinates.


function trajOut = trajMirror(trajData, coordCols)

if isa(trajData,'double')
    trajData = {trajData};
end

if any(cellfun(@isempty, trajData))
    error('Cell array must not have empty cells')
end

if nargin < 2
    coordCols = [1 2];
end

tmpDat = cellfun(@(tr)  tr(:,coordCols) ,trajData,'uniformoutput',0);

% Rotate onto axis 2 (y) and shift to origin.
[tmpDat, trajStart, reverseRefDirs]  = trajRot(tmpDat,2,[1 2],0,1);

% Flip
tmpDat = cellfun(@(tr) [tr(:,1)*-1,tr(:,2)],tmpDat,'uniformoutput',0);

% Reverse rotation and translation
tmpDat = trajRevRot(tmpDat, 2,[1 2],reverseRefDirs,trajStart);

for i = 1:numel(tmpDat)       
    
    trajData{i}(:,coordCols(1)) = tmpDat{i}(:,1);
    trajData{i}(:,coordCols(2)) = tmpDat{i}(:,2);
    
end

trajOut = trajData;



