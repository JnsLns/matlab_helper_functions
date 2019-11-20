
% function [trajOut, trajStart, reverseRefDirs] = trajRot(trajData, rotateToAxis, coordCols, refDirs, shiftToOrigin)
%
% Rotates trajectories of one- to three-dimensional data points such that
% after rotation a reference axis defined relative to the trajectory is
% parallel to one of the coordinate axes. The reference axis can be speci-
% fied freely for each trajectory but is by default defined as the vector
% between the first and last data point of the trajectory, so that in effect
% the trajectory is rotated such that it starts and ends on a line parallel
% to the target coordinate axis. If a matrix in trajData contains only one
% data point, the coordinate origin is assumed as first data point, and the
% input point as the last. (note: Arbitrary rotation angles can be
% achieved through defining custom refDirs; Reversing the entire transformation
% including shift and rotation can be realized by feeding the outputs
% trajStart  (holding the original starting points of the trajectories) and
% reverseRefDirs (holding the ideal path of the rotated trajectories) to the
% function trajRotRev, along with the same value of rotateToAxis as used
% in the call to trajRot to trajRevRot. Note that this reverse rotation can
% be applied to any data point, e.g. allowing to rectify trajectories, com-
% pute data points using the rectified versions, and then transfrom the new-
% ly derived points back into the original coordinate space, for further
% analysis, plotting etc.).
% IMPORTANT: Be aware that the default rotation of this function is
% left-handed and such that the last trajectory point ends up on the
% *positive* y-axis - which makes sense only if trajectories in the initial
% coordinate frame travel into positive y-direction as well (whereas if
% trajectories travel from positive to negative y-values, this effectively
% results in the trajectory being mirrored across the reference axis!). 
%
% trajData: n-by-m cell array holding in each element one trajectory (max.
% three dimensions) in the form of an n-by-m matrix: In this matrix,
% n = data points, and m = point coordinates. trajData may also be supplied
% as an n-by-m matrix of that form holding only one trajectory.
%
% rotateToAxis (optional, default 1): number of axis to which reference
% vector (refDirs) should be parallel after rotation; the number refers to
% the order of dimensions given in coordCols (e.g., if x-coordinates are in
% column one, and y-coordinates in column two (coordCols = [1 2]) and
% trajectories should be rotated to the y-axis, rotateToAxis should be 2).
%
% coordCols (optional, default [1 2 3]): n-element row vector giving the
% column numbers of those columns in each trajectory matrix that hold
% coordinates (e.g., x, y, z).
%
% refDirs (optional, default is to compute refDir individually for each
% trajectory as the vector between its first and last data point):
% n-by-m vector holding in each row a vector that defines a reference
% direction according to which the rotation angle is determined.
% In effect, trajectory_i is rotated by the angle between refDirs_i and the
% coordinate axis denoted by rotateToAxis. Column number of refDirs, m,
% should be equal to numel(coordCols). Row number, n, of refDirs should be
% numel(trajData) or one. If size(refDirs,1) == 1, the reference
% direction it holds is used for all trajectories. The same is true if instead
% 1 < n ~= numel(trajData). 
% refDirs == [0] overrides all other behaviors and is equal to default.
%
% shiftToOrigin (optional, default 0): If 1, the trajectory is not only
% rotated, but it is is also translated such that its first data point lies
% in the origin of the coordinate frame. 

function [trajOut, trajStart, reverseRefDirs] = trajRot(trajData, rotateToAxis, coordCols, refDirs, shiftToOrigin)

if isa(trajData,'double')
    trajData = {trajData};
end

if any(cellfun(@isempty, trajData))
    error('Cell array must not have empty cells')
end

if nargin < 2
    rotateToAxis = 1;
end

if nargin < 3
    coordCols = [1 2 3];
end

if nargin < 5
    shiftToOrigin = 0;
end

origin =  zeros(1,numel(coordCols));

% ensure that trajData has only one column (avoids problems with cellfun
% outputs, permute, etc.)
trajData = reshape(trajData,numel(trajData),1);

% Check whether any trajectories consist of only one data point. If so, add
% origin as first data point (i.e., row of zeros at the top). These dummy
% points are removed later.
isSingle = cellfun(@(x) size(x,1)==1,trajData);
trajData = cellfun(@(x) [zeros(1*(size(x,1)==1),size(x,2));x],trajData,'UniformOutput', 0);

% Get first data point of each trajectory as output.
trajStart = cellfun(@(x) x(1,coordCols), trajData, 'UniformOutput', 0);

% Store original coordinates of first data point to later reverse shift
trajStartMem = cellfun(@(x) x(1,coordCols), trajData, 'UniformOutput', 0);
 
% Shift all trajectories such that the first point lies at the origin
% of the coordinate system.
trajData = trajShift(trajData, 1, origin, coordCols);

% unit vector along the target axis
unitVec = [0 0 0];
unitVec(rotateToAxis) = 1;

% how many dimensions less than three does the data have? (used for padding
% to three dims for computations later on)
missingDims = 3-numel(coordCols);

% loop over trajectories and rotate each
for curTraj = 1:numel(trajData)
    
    % Get one trajectory from trajData to work on
    traj = trajData{curTraj}(:,coordCols);
    
    % Determine reference direction vector that is used to determine rotation angle.
    if nargin < 4 || (numel(refDirs) == 1  && refDirs == 0)        
        
        % If no other reference vector is supplied, use the general direction of
        % the current trajectory, i.e., the vector between last and first point.
        refDir = traj(end,:) - traj(1,:);
                
    elseif nargin >= 4       
        
        % If refDirs are supplied, get the one for this trajectory.
        % If supplied but different number than trajectories, use the first
        % one for all trajectories.
        if size(refDirs,1) == numel(trajData)
            refDir = refDirs(curTraj,:);
        else            
            refDir = refDirs(1,:);
            warning('in function trajRot: number of supplied reference directions (refDirs) is different from number of trajectories. Using first row of refDirs.');
        end                          
        
    end             
             
    % if data has less than three dimensions (columns), pad to three dims
    % with zero columns; same for refDir
    traj(:,end+1:end+missingDims) = 0;
    refDir(:,end+1:end+missingDims) = 0;
            
    % compute rotation axis rotax as (unit) normal of the plane spanned by
    % target coordinate axis (unit vector) and the to-be-rotated vector refDir   
    rotax = cross(refDir,unitVec);
    if norm(rotax) ~= 0 % to prevent Nans through division by zero
        rotax = rotax/norm(rotax); 
    end    
    rotaxes{curTraj,1} = rotax;
    
    % compute theta as (left-handed) angle between target axis and
    % reference direction (and store for output)
    theta  = acos(dot(refDir,unitVec)/(norm(refDir)));
    thetas{curTraj,1} = theta;
    
    % Rotate using Rodrigues' rotation formula (v is the current data point
    % being rotated). The anonymous function implements the rotation formula
    % and is applied to each row (data point) of the current trajectory.
    rodrot = @(v) (v*cos(theta) + (cross(rotax,v)) * sin(theta) + rotax * (dot(rotax,v)) * (1 - cos(theta)));
    traj = cell2mat(arrayfun(@(trRows) rodrot(traj(trRows,:)), 1:size(traj,1), 'UniformOutput', 0)');
    
    % replace data in trajData (removing padded dimensions in the process)
    trajData{curTraj}(:,coordCols) = traj(:,1:end-missingDims);
    
end

if exist('thetas','var')
    % use the stored thetas (rotation angles) to further rotate each trajectories'
    % refDir starting with its *already rotated* form, to obtain reverseRefDirs
    % that can be used to reverse the rotation performed by the current function when
    % employed in a call to trajRevRot. Note that each refDir, rotated by theta,
    % is parallel to the target axis rotateToAxis so that we can simply use
    % the unit vector unitVec and rotate it again by theta to get reverseRefDir.
              
    if nargout >= 3                            
        reverseRefDirs = repmat(unitVec,numel(trajData),1);        
        reverseRefDirs = mat2cell(reverseRefDirs,ones(size(reverseRefDirs,1),1));       
        rrot = @(v,the,rax) (v*cos(the) + (cross(rax,v)) * sin(the) + rax * (dot(rax,v)) * (1 - cos(the)));
        reverseRefDirs = cell2mat(cellfun(@(rrds,ts,ras) rrot(rrds,ts,ras), ...
            reverseRefDirs, thetas, rotaxes, 'UniformOutput', 0));
        reverseRefDirs = reverseRefDirs(:,1:numel(coordCols));
    end
else
    reverseRefDirs = [];
end

% Shift all trajectories' starting points back to their original position
if ~shiftToOrigin
    trajData = trajShift(trajData, 1, cell2mat(trajStartMem), coordCols);
end

% remove dummy rows in "single-point trajectories" and return
trajOut = cellfun(@(tr,is)  tr(is+1:end,:),trajData,num2cell(isSingle),'uniformoutput',0);

end

