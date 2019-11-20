% 
% function angleToIorCell = trajAngleToIorMinTgtMinConfound_WIP(trajData,iorPos,confoundPos,tgtPos,coordCols,outCol)
% %
% % This function is supposed to compute the angular deviation at each time step
% % of a trajectory from the direction toward a given item position (iorPos) 
% %
% % HEAVILY WIP!!!!! NOT WORKING PROPERLY!!! HAVE TO COMPUTE *DIRECTED*
% % ANGLE !!!
% %
% % THE REMAINING DOC IS OUTDATED. THIS FUNCTION HAS BEEN CHANGED!!!
% % 
% %
% % This function computes the magnitude of the velocity component in the
% % direction of a given reference point, for each timestep of a 2d-trajectory
% % or for multiple trajectories in a cell array. In effect, the original
% % coordinate system of the input data is rotated such that the y-axis is
% % parallel to a straight line between the first point of a given input
% % velocity vector and a given reference point, while the input velocity
% % vector stays static, and the y-component of the old vector in the new
% % coordinate system is returned. The reference line determining coordinate 
% % rotation is computed anew for each timestep.
% %
% % Depending on the value of outCol, the output of the function is either 
% % a cell array identical to the input array but with the velocity
% % information added as a column to each trajectory matrix, or a cell array
% % with each cell holding a column vector with only the velocity data of the
% % respective trajectory.
% %
% % Note: The velocity value in a given row of the output data refers to the
% % trajectory segment that starts at the time point given in the same row and
% % ending at the time point given in the subsequent row. The last row of the
% % output data column is set to 0.
% %
% % 
% % trajData: n-by-m cell array holding in each element one trajectory in the
% % form of an n-by-m matrix: In this matrix, n = data points, and m = point
% % coordinates. trajData may also be supplied as an n-by-m matrix of that
% % form holding only one trajectory. 
% %
% % iorPos: n-by-2 matrix where n = numel(trajData), each row providing a
% % two-coordinate reference point for the trajectory at the corresponding 
% % linear index in trajData. For each timestep of that trajectory, the
% % velocity toward that point will be computed. Regardless of numel(trajData),
% % iorPos may alternatively provided as a 2-element row vector, in which
% % case that reference point is used for all trajectories.
% %
% % coordCols (optional, default [1 2]): 2-element row vector giving the
% % column numbers of those columns in the trajectory matrices that hold
% % position coordinates.
% % 
% % tCol (optional, default 3): integer giving the column number of the
% % column in the trajectory matrices that holds time values. If set to 0,
% % the movement distance in each time step is not normalized by time, so that
% % the function returns the raw movement distance in the direction of the
% % reference position instead of the velocity.
% %
% % outCol (optional, default 'add'): Determines where the output column of
% % velocity values goes. If 'add', it is added as an additional column
% % to the input trajectory matrices. If 'single', only the resulting column is
% % returned, without the original input data. Alternatively, outCol may be
% % a scalar value specifying a column in the input trajectory matrices, which
% % is replaced by the computed values.
% %
%   
% %msgbox('THIS FUNCTION IS WIP AND DOES NOT WORK PROPERLY!')
% warning('trajangleToIor is WIP / PROTOTYPE and might not work as intended!!!')
% 
% if isa(trajData,'double')
%     trajData = {trajData};
% end
% 
% if any(cellfun(@isempty, trajData))
%     error('Cell array must not have empty cells')
% end
% 
% if nargin < 5
%     coordCols = [1 2];
% end
% 
% if size(iorPos,1) < numel(trajData)
%     iorPos = repmat(iorPos,numel(trajData),1);
%     confoundPos = repmat(confoundPos,numel(trajData),1);
%     tgtPos = repmat(tgtPos,numel(trajData),1);
% end
%     
% 
% for curTraj = 1:numel(trajData)
% 
%     tr = trajData{curTraj};           
%         
%     angleToIor = zeros(size(tr,1)-1,1);         
%     angleToConfound = zeros(size(tr,1)-1,1);         
%     angleToTgt = zeros(size(tr,1)-1,1);         
%     
%     if nargin < 6  || strcmp(outCol,'add')
%         outCol = size(tr,2)+1;
%     end
%     
    for curRow = 1:size(tr,1)-1
       
        % vectors toward items computed from each trajectory data points
        % movement vector
        mv = tr(curRow+1,coordCols) - tr(curRow,coordCols);
        % item of interest vector 
        iv = iorPos(curTraj,:) - tr(curRow,coordCols);        
        % confounding item
        cv = confoundPos(curTraj,:) - tr(curRow,coordCols);        
        % target item (used to compute direct path from each data point)
        tv = tgtPos(curTraj,:) - tr(curRow,coordCols);        
                
        % left-handed angle to item of interest       
        angleToIor(curRow,1) = deg2rad(atan2d(iv(1)*mv(2)-iv(2)*mv(1),iv(1)*mv(1)+iv(2)*mv(2)));                                              
        
        % left-handed angle to confound
        angleToConfound(curRow,1) = deg2rad(atan2d(cv(1)*mv(2)-cv(2)*mv(1),cv(1)*mv(1)+cv(2)*mv(2)));
        
        % left-handed angle to target
        angleToTgt(curRow,1) = deg2rad(atan2d(tv(1)*mv(2)-tv(2)*mv(1),tv(1)*mv(1)+tv(2)*mv(2)));
                        
    end
%     
%     outAngle = angleToIor - angleToConfound
%     
%     dummy = 0;
%     if strcmp(outCol,'single')
%         trajData{curTraj} = [angleToIor-angleToTgt;dummy];
%     else        
%         trajData{curTraj}(:,outCol) = [angleToIor-angleToTgt;dummy];
%     end
%                     
% end
% 
% angleToIorCell = trajData;
% 
