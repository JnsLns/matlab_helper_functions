
function speedCell = trajSpeed(trajData,coordCols,tCol,outCol)
%
% This function computes the speed for each timestep of a 2d-trajectory or
% for multiple trajectories provided as matrices in a cell array.
% Depending on the value of outCol, the output of the function is either 
% a cell array identical to the input array but with the speed
% information added as a column to each trajectory matrix, or a cell array
% with each cell holding a column vector with only the speed data of the
% respective trajectory.
%
% Note: The speed value in a given row of the output data refers to the
% trajectory segment that starts at the time point given in the same row and
% ending at the time point given in the subsequent row. The output matrices 
% are one row shorter than the input matrices, missing the last row.
%
% trajData: n-by-m cell array holding in each element one trajectory in the
% form of an n-by-m matrix: In this matrix, n = data points, and m = point
% coordinates. trajData may also be supplied as an n-by-m matrix of that
% form holding only one trajectory. 
%
% coordCols (optional, default [1 2]): 2-element row vector giving the
% column numbers of those columns in the trajectory matrices that hold
% position coordinates.
% 
% tCol (optional, default 3): scalar giving the column number of the
% column in the trajectory matrices that holds time values. If set to 0,
% the movement distance in each time step is *not* normalized by time, so
% that the function returns movement distance for each timestep instead of
% speed.
%
% outCol (optional, default 'add'): Determines where the output column of
% speed values goes. If 'add', it is added as an additional (trailing) column
% to the input trajectory matrices. If 'single', only the resulting column is
% returned, without the original input data. Alternatively, outCol may be
% a scalar value specifying a column in the input trajectory matrices, which
% is replaced by the computed values.
%
  
if isa(trajData,'double')
    trajData = {trajData};
end

if any(cellfun(@isempty, trajData))
    error('Cell array must not have empty cells')
end

if nargin < 2
    coordCols = [1 2];
end

if nargin < 3
    tCol = 3;
end
  

for curTraj = 1:numel(trajData)

    tr = trajData{curTraj};           
        
    speed = zeros(size(tr,1)-1,1);         
    
    if nargin < 4  || strcmp(outCol,'add')
        outCol = size(tr,2)+1;
    end
    
    for curRow = 1:size(tr,1)-1                       
        
        % Get movement distance
        movementDistance = norm(tr(curRow+1,coordCols) - tr(curRow,coordCols));
   
        % Get time interval
        if tCol ~= 0
            deltat = (tr(curRow+1,tCol)-tr(curRow,tCol));
        else
            deltat = 1;
        end
        
        speed(curRow) = movementDistance/deltat;                                              
        
    end
        
    if strcmp(outCol,'single')
        trajData{curTraj} = speed;
    else
        trajData{curTraj}(end,:) = []; % remove final row from original data
        trajData{curTraj}(:,outCol) = speed;
    end
                    
end

speedCell = trajData;

