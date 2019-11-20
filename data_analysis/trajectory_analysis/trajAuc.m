% function [AUCs] = trajAuc(trajData, xyCols, absFlag, onlyLR)
% 
% Computes the area under curve (AUC) for a set of 2d trajectories (or other
% 2d functions), each provided as an individual matrix in a cell array.
% It returns AUCs, a column vector with each row holding the AUC value for
% the trajectory at the corresponding linear index in the cell array.
% AUC is computed as integral of the trajectory over the straight line
% connecting its first and last data point.
% Note that deviations to the left of this direct path are counted 
% as negative and deviations to the right are counted as positive.
% Thus, for instance, a positive AUC only indicates *more* or overall
% farther deviation to the right side, and not necessarily deviation
% exclusively to the right. Similarly, zero AUC does not necessarily indi-
% cate no deviation from the direct path, but may as well result from
% equally far deviations to either side. Setting onlyLR allows to count only 
% one side of deviation. Absolute AUCs (counting both left and right
% deviations as positive) may be obtained by setting absFlag to 1.
% 
% ***Note that valid results can only be expected when the trajectory data
% points are defined over strictly monotonically increasing sample points.***
% This means it is mostly not a good idea to use this function for real
% space curves that may turn backwards at some point; ideally use it for
% data points defined over time instead (i.e., supply monotonically increa-
% sing time points as in the column denoted by xyCols(2)).
%
%
% ----IMPORTANT NOTE ABOUT ABSFLAG AND onlyLR---
% Absolute AUCs (absFlag == 1) are computed by simply making all data
% points positive before numerical integration. This procedure may lead to
% imprecisions in total area if there are no data points at or near axis
% crossings. These are usually small but may become larger when
% trajectories consist of a small number of data points (in that case,
% trajectories should be interpolated before subjecting them to this
% function). The same is true for onlyLR (if 1 or 2, negative or positive
% data points, respectively, are simply set to zero before integrating).
% 
%
% trajData: n-by-m cell array holding in each element one trajectory in the
% form of an p-by-q matrix: In this matrix, p = data points, and q = point
% coordinates. trajData may also be supplied as an p-by-q matrix of that
% form holding only one trajectory.
%
% xyCols (optional, default [1 2]): 2-element row vector giving the
% column numbers of those columns in each trajectory matrix that hold
% the x and y coordinates of the trajectory. Note that, for the mapping of 
% negative AUCs to leftward deviation and vice versa to be true, this
% mapping must be consistent with that used in other computations and when
% trajectories were acquired.
%
% absFlag (optional, default 0): If 0, AUCs are computed with data points
% to the right of the direct path between first and last data point
% counting as positive and vice versa. If 1, the AUC is computed as the sum
% of the absolute distance values (so that the AUC gives a measure of
% overall deviation from the direct path, disregarding in which direction).
% 
% onlyLR (optional, default 0, i.e., deactivated): If 1, only deviations to
% the left of the direct path are counted (as positive in this case). If 2,
% only deviations to the right of the direct path are counted (as well as
% positive). Read note above!


function [AUCs] = trajAuc(trajData, xyCols, absFlag, onlyLR)

if nargin < 4
    onlyLR = 0;
end

if nargin < 3
    absFlag = 0;
end

if nargin < 2
    xyCols = [1 2];
end

if isa(trajData,'double')            
    trajData = {trajData};            
end

if any(cellfun(@isempty, trajData))
    error('Cell array must not have empty cells')
end

% ensure that trajData has only one column (avoids problems with cellfun
% outputs, permute, and in case of interaction with other functions)
trajData = reshape(trajData,numel(trajData),1);

% Check that second coordinate column always is monotonically increasing
if ~all(cellfun(@(tr) all(diff(tr(:,xyCols(2)))>=0),trajData))
    error('Data in column xyCols(2) must be monotonically increasing for all trajectories.')
end

% Make parallel to y-axis and shift to origin
trajData = trajRot(trajData,2,xyCols,0,1);

% count only left- or rightward deviation if desired
if onlyLR == 1 % count only left (as positive)    
    for curTr = 1:numel(trajData)        
        trajData{curTr}(trajData{curTr}(:,xyCols(1)) > 0,1) = 0;        
    end
    absFlag = 1;
elseif onlyLR == 2 % count only right (as positive)    
    for curTr = 1:numel(trajData)        
        trajData{curTr}(trajData{curTr}(:,xyCols(1)) < 0,1) = 0;        
    end    
    absFlag = 1;
end
    
% compute area using trapz
if ~absFlag
    AUCs = cellfun(@(x) trapz(x(:,xyCols(2)),x(:,xyCols(1))),trajData);
elseif absFlag == 1
    AUCs = cellfun(@(x) trapz(x(:,xyCols(2)),abs(x(:,xyCols(1)))),trajData);
end


