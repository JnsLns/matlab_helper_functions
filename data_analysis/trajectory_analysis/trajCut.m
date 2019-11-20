% function trajOut = trajCut(trajData,fromToRows,fromToTimes,tCol)
%
% Cuts off the first and last portions of all trajectories in cell array
% trajData, according to start and end points defined by row numbers in
% fromToRows or by first and last time points in fromToTimes. Returns a
% cell array similar to trajData, containing the cut trajectories instead.
%
% trajData: n-by-m cell array holding in each element one trajectory in the
% form of an n-by-m matrix: In this matrix, n = timepoints, and m = types
% of data (x,y,t,...; there may be an arbitrary number of columns). trajdata
% may also be supplied as an n-by-m matrix of that form, holding only one
% trajectory.
%
% fromToRows: n-by-2 matrix, where n is equal to the number of elements in
% trajData, with each row holding the row numbers of the desired first and
% last point for the trajectory with the corresponding linear index in
% trajData. All rows *before and after* these points are removed
% (*excluding* first and last rows themselves). If second element of fromToRows
% is 0, then only the trajectory portion up to fromToRows's first element
% is removed. fromToRows is completely disregarded if fromToTimes is defined.
%
% fromToTimes (optional, overrides fromToRows if defined): n-by-2 matrix,
% where n is equal to the number of elements in trajData, with each row
% holding the time values of the desired first and last point for the
% trajectory with the corresponding linear index in trajData. 
%
% tCol (onyl needed if fromToTimes is defined, default 3): column number of
% time values in matrices in trajData.



function trajOut = trajCut(trajData,fromToRows,fromToTimes,tCol)

if isa(trajData,'double')
    trajData = {trajData};
end

if any(cellfun(@isempty, trajData))
    error('Cell array must not have empty cells')
end

% ensure that trajData has only one column (avoids problems with cellfun
% outputs, permute, and in case of interaction with other functions)
trajData = reshape(trajData,numel(trajData),1);

trajOut = cell(size(trajData));

% find rows in trajectory matrices that correspond to time values provided
% in fromToTimes
if nargin > 2   
    if nargin < 4
        tCol = 3;
    end    
    fromToRows = cell2mat(cellfun(@(tr,tFirst,tLast)  [find(tr(:,tCol)==tFirst),find(tr(:,tCol)==tLast)], ...
    trajData,num2cell(fromToTimes(:,1)),num2cell(fromToTimes(:,2)),'uniformoutput',0));           
end

% go through trajectories and cut off start and end rows
for currentTraj = 1:numel(trajData)

    traj = trajData{currentTraj};    
    from = fromToRows(currentTraj,1);        
        if fromToRows(currentTraj,2) ~= 0
        to = fromToRows(currentTraj,2);
    else
        to = size(traj,1);
    end    
    trajOut{currentTraj} = traj(from:to,:);
        
end

end