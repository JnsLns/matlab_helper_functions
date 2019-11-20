
function velTowPointCell = trajVelTowPoint(trajData,refPos,coordCols,tCol,outCol)
%
% This function computes the magnitude of the velocity component in the
% direction of a given reference point, for each timestep of a 2d-trajectory
% or for multiple trajectories in a cell array. In effect, the original
% coordinate system of the input data is rotated such that the y-axis is
% parallel to a straight line between the first point of a given input
% velocity vector and a given reference point, while the input velocity
% vector stays static, and the y-component of the old vector in the new
% coordinate system is returned. The reference line determining coordinate 
% rotation is computed anew for each timestep.
%
% Depending on the value of outCol, the output of the function is either 
% a cell array identical to the input array but with the velocity
% information added as a column to each trajectory matrix, or a cell array
% with each cell holding a column vector with only the velocity data of the
% respective trajectory.
%
% Note: The velocity value in a given row of the output data refers to the
% trajectory segment that starts at the time point given in the same row and
% ending at the time point given in the subsequent row. The last row of the
% output data column is set to 0.
%
% 
% trajData: n-by-m cell array holding in each element one trajectory in the
% form of an n-by-m matrix: In this matrix, n = data points, and m = point
% coordinates. trajData may also be supplied as an n-by-m matrix of that
% form holding only one trajectory. 
%
% refPos: n-by-2 matrix where n = numel(trajData), each row providing a
% two-coordinate reference point for the trajectory at the corresponding 
% linear index in trajData. For each timestep of that trajectory, the
% velocity toward that point will be computed. Regardless of numel(trajData),
% refPos may alternatively provided as a 2-element row vector, in which
% case that reference point is used for all trajectories.
%
% coordCols (optional, default [1 2]): 2-element row vector giving the
% column numbers of those columns in the trajectory matrices that hold
% position coordinates.
% 
% tCol (optional, default 3): integer giving the column number of the
% column in the trajectory matrices that holds time values. If set to 0,
% the movement distance in each time step is not normalized by time, so that
% the function returns the raw movement distance in the direction of the
% reference position instead of the velocity.
%
% outCol (optional, default 'add'): Determines where the output column of
% velocity values goes. If 'add', it is added as an additional column
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

if nargin < 3
    coordCols = [1 2];
end

if nargin < 4
    tCol = 3;
end

if size(refPos,1) < numel(trajData)
    refPos = repmat(refPos,numel(trajData),1);
end
    

for curTraj = 1:numel(trajData)

    tr = trajData{curTraj};           
        
    velTowPoint = zeros(size(tr,1)-1,1);         
    
    if nargin < 5  || strcmp(outCol,'add')
        outCol = size(tr,2)+1;
    end
    
    for curRow = 1:size(tr,1)-1
       
        movementVector = tr(curRow+1,coordCols) - tr(curRow,coordCols);
        
        % Reference direction for rotation is defined by line between
        % target object and first point defining current movement vector.
        referenceVector = refPos(curTraj,:) - tr(curRow,coordCols);        
        rotatedMovement = cell2mat(trajRot(movementVector,2,1:numel(coordCols),referenceVector,1));
        
        % use y-value of rotated movement, which equals movement toward the
        % target object in the current timestep, and divide by length of
        % the timestep to obtain velocity.        
        if tCol ~= 0
            deltat = (tr(curRow+1,tCol)-tr(curRow,tCol));
        else
            deltat = 1;
        end
        velTowPoint(curRow,1) = rotatedMovement(2) ./ deltat;                           
        
    end
    
    dummy = 0;
    if strcmp(outCol,'single')
        trajData{curTraj} = [velTowPoint;dummy];
    else        
        trajData{curTraj}(:,outCol) = [velTowPoint;dummy];
    end
                    
end

velTowPointCell = trajData;

