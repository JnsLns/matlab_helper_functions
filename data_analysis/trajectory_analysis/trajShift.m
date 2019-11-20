% function trajOut = trajShift(trajData,referenceSample,shiftToPoint,coordCols)
%
% Translates an n-dimensional trajectory according to the vector between a
% reference point that is on the trajectory, and a target point supplied
% freely. Returns a cell array similar to trajData i.e., containing
% trajectory matrices with the same colum mapping but with transformed
% coordinate data.
%
% trajData: n-by-m cell array holding in each element one trajectory in the
% form of an n-by-m matrix: In this matrix, n = data points, and m = point
% coordinates. trajData may also be supplied as an n-by-m matrix of that
% form holding only one trajectory.
%
% referenceSample (optional): n-element column vector, for each trajectory
% giving the number of the row from which the reference data should be taken
% (i.e., coordinates of the point that will end up on the target point). If
% referenceSample has less elements than trajData (e.g., one), then its
% first element is used for all trajectories. Default is [1].
%
% shiftToPoint (optional): n-by-m matrix where n is the number of
% trajectories and m is the number of dimensions. Each row holds the point
% coordinates to which the trajectory in the corresponding element of
% trajData will be shifted. If n is not equal to the number of
% trajectories, then the first row is used for all trajectories. Default is
% [0 0 0].
%
% coordCols (optional): n-element row vector giving the column numbers of
% those columns in each trajectory matrix that hold coordinates (e.g., x,
% y,z). Default is [1 2 3].


function trajOut = trajShift(trajData,referenceSample,shiftToPoint,coordCols)


if isa(trajData,'double')
    trajData = {trajData};
end

if any(cellfun(@isempty, trajData))
    error('Cell array must not have empty cells')
end

if nargin < 4
    coordCols = [1 2 3];
end

if nargin < 3
    shiftToPoint = [0 0 0];
end

if nargin < 2
    referenceSample = 1;    
end

% ensure that trajData has only one column (avoids problems with cellfun
% outputs, permute, and in case of interaction with other functions)
trajData = reshape(trajData,numel(trajData),1);

% there should be as many referenceSample elements as trajectories.
% If not, use first element for all trajectories.
if numel(referenceSample) ~= numel(trajData)
    referenceSample = referenceSample(1);
end
if numel(referenceSample) == 1
    referenceSample = repmat(referenceSample,numel(trajData),1);
end

% if there is a different number of shiftToPoint's than trajectories, use
% first row of shiftToPoint for all trajectories.
if size(shiftToPoint,1) ~= numel(trajData)
    shiftToPoint = shiftToPoint(1,:);
end
if size(shiftToPoint,1) == 1
    shiftToPoint = repmat(shiftToPoint,numel(trajData),1);
end

% prepare output cell array
trajOut = trajData;

% go through trajectories and shift
for currentTraj = 1:numel(trajData)
    
    traj = trajData{currentTraj}(:,coordCols);
    
    % Get reference point on trajectory
    refPoint = traj(referenceSample(currentTraj),:);
    % Get target point
    tgtPoint = shiftToPoint(currentTraj,:);
    % Compute shift vector
    shiftVec = tgtPoint - refPoint;
    % Apply shift
    traj = bsxfun(@plus,traj,shiftVec);
    % Store in output cell
    trajOut{currentTraj}(:,coordCols) = traj;
    
end


end