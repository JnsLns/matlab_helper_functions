function [k,interCurve,arcSpanned] = curvatureOsc(vertices,interpToUniformSegmentLengthFlag,interpInterval,symThresh,symTimeThresh,roundToDec)
% [k,interCurve,arcSpanned] = curvatureOsc(vertices,interpToUniformSegmentLengthFlag,interpInterval,symThresh,symTimeThresh,roundToDec)
%
% Returns curvature of an ordered list of vertices that define a plane curve
% (m-by-2 matrix where m are vertices and columns are x,y coordinates).
% The returned measure is the reciprocal radius of the osculating circle
% through a vertex and its neighboring two vertices (nan is returned for
% first and last data point). Zero angle at a vertice yields a curvature 
% of inf. 
%
% interpToUniformSegmentLengthFlag (optional, default 0): If set to 1,
% before curvature for each vertex is assessed, the curve is linearly
% interpolated such that each segment of the resulting curve has the same
% length. That way, the impact of segment length (i.e., vertex distances)
% is removed from the curvature measure (With the osculating circle method,
% the same angle between two successive curve segments may results in vastly
% different curvature values, depending on the segment lengths as a major
% factor impacting the distance of the current to the neighboring two
% vertices and thus circle size; whether it makes sense to remove that
% effect depends on the type of data and analysis at hand).
% If enabled, interpInterval must as well be specified, as a scalar that
% determines the distance between vertices in the interpolated curve. The
% returned vector of vertex curvatures k will have a length equal to the
% number of points in the resulting interpolated curve. Note that the
% segment length will have major impact on curvature and should be adjusted
% such that the returned curvature actually reflects the aspects of the
% considered curves that is relevant in the context at hand. While
% curvature for a vertex on a straight line of neighboring points will
% always be zero (regardless of the points distances), the highest
% achievable curvature value is dependent on the segments length: For
% instance, if segment length is 5, then the minimum possible radius of an
% osculating circle using three points is 2.5, which equals a highest 
% achievable curvature value of of 1/2.5 = 0.4.
%
% The remaining arguments are passed to interpToUniformSegments, which is
% for interpolation if the flag is set; see INTERPTOUNIFORMSEGMENTS for
% details. They are optional arguments with default values (here):
% symthresh: 0.1
% symTimeThresh: 0.25 seconds
% roundToDec: -1 (i.e., rounding disabled, increase>1 for performance at
%                 cost of precision)
%
% Optional output interCurve is the interpolated version of the curve, in
% same format as input vertices.
%
% Optional output arcSpanned gives the arc length of the original curve
% that is spanned by (i.e., lies in between the two vertices of) each segment
% of the interpolated curve.
%
% See also INTERPTOUNIFORMSEGMENTS.

if nargin < 2
    interpToUniformSegmentLengthFlag = 0;
end

if nargin == 2 && interpToUniformSegmentLengthFlag
    error('If interpToUniformSegmentLengthFlag is set to 1, interpInterval must be specified as well.');
end

if nargin < 4
    symThresh = 0.1;
end

if nargin < 5
    symTimeThresh = 0.25;
end

if nargin < 6
    roundToDec = -1;
end


if interpToUniformSegmentLengthFlag
    [interCurve,arcSpanned] = interpToUniformSegments(vertices(:,1),vertices(:,2),interpInterval,0,symThresh,symTimeThresh,roundToDec);            
    vertices = interCurve;
end

% Prepare curvature vector
k = zeros(size(vertices,1)-2,1);
k = [nan;k;nan];

for curRow = 2:size(vertices,1)-1
    
    % Get current vertices and vertices around it
    x1 = vertices(curRow-1,1);
    y1 = vertices(curRow-1,2);
    x2 = vertices(curRow,1);
    y2 = vertices(curRow,2);
    x3 = vertices(curRow+1,1);
    y3 = vertices(curRow+1,2);
            
    % Compute circle through these points
    % mid point y-coordinate
    y_m = ((x3^2 - x1^2 + y3^2 - y1^2) * (x2-x1) - (x2^2 - x1^2 + y2^2 - y1^2) * (x3-x1)) / ...
        (2*((y3-y1)*(x2-x1) - (y2-y1)*(x3-x1)));
    % mid point x-coordinate        
    x_m = ((x2^2-x1^2) + (y2^2-y1^2) - 2*y_m*(y2-y1)) / (2*(x2-x1));
    if isnan(x_m) || isinf(x_m)
    x_m = ((x3^2-x1^2) + (y3^2-y1^2) - 2*y_m*(y3-y1)) / (2*(x3-x1));
    end    
    % radius
    r = sqrt((x1-x_m)^2 + (y1-y_m)^2);
        
    % curvature is reciprocal of radius
    if ~isnan(1/r ) 
        k(curRow) = 1/r;
    elseif sign(dot([x2,y2]-[x1,y1],[x3,y3]-[x2,y2])) == 1 % edges are parallel
        k(curRow) = 0;
    elseif sign(dot([x2,y2]-[x1,y1],[x3,y3]-[x2,y2])) == -1 % edges are anti-parallel
        k(curRow) = inf;
    end
    
    
end
    


