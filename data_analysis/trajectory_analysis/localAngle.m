function [maxAngle, localAngles] = localAngle(x, y, segLen, stepSize, symThresh, symTime, decRound)
% function [maxAngle, localAngles] = localAngle(x, y, segLen, stepSize, symThresh, symTime, decRound)
%
% Compute local angles along a trajectory according to the algorithm
% described in Lins & Schöner (2019, Atten Percept Psychophys, Appendix B).
%
% The algorithm works by fitting the three vertices of a two-segment
% line fragment with segment length segLen to successive points along the
% trajectory, separated by steps of stepSize arc length (both parameters
% in units of x, y) and computing the resulting angle (in radians) between
% the fragment's first and second segment at each trajectory position.
% Obtained angles are between 0 radians (straight line) and pi radians
% (antiparallel segments). Parameter stepSize essentially determines the
% sampling density whereas segLen influences the scale of trajectory turns
% that will or will not be registered as an angled portion (e.g., large
% segLen will lead to smaller sharp turns being disregarded). These
% parameters should be adjusted to the nature of the data at hand and to
% the aim of the procedure.
%
% Note: For the paper mentioned above parameters were segLen = 15, stepSize
%       = 1, symThresh = 0.05, symTime = 1, decRound = 2.
%
% __Input arguments__
%
% x,y           Trajectory coordinate vectors, see above.
%
% segLen        Segment length of fitted fragment. (internally passed to
%               interpToUniformSegments_fixedn)
%
% stepSize      Steps of trajectory arc length at which fragment is fitted.
%               (internally passed to interpToUniformSegments_fixedn)
% 
% symThresh     Optional, default 0.05.
%               Spatial threshold to switch to symbolic math mode within
%               interpToUniformSegments_fixedn. (internally passed to
%               interpToUniformSegments_fixedn)
%
% symTime       Optional, default 1 second.
%               Time threshold to switch to symbolic math mode within
%               interpToUniformSegments_fixedn. (internally passed to
%               interpToUniformSegments_fixedn)
%
% decRound      Optional, default 2.
%               Determines to how many decimals coordinate input is rounded
%               prior to computing interpolation in
%               interpToUniformSegments_fixedn. (internally passed to
%               interpToUniformSegments_fixedn) 
%   
% __Output arguments__
%
% maxAngle      Maximum angle obtained anywhere along the trajectory. 
%
% localAngles   Vector containing all angles obtained along the trajectory.
%               localAngles(1) is the angle obtained when placing the first 
%               vertex of the fitted fragment at the trajectory's first
%               data point. Subsequent elements in localAngles are the
%               angles obtained for placing it at successive positions
%               separated by stepSize steps of trajectory arc length. Note
%               that computation stops as soon as the third vertex of the
%               fragment cannot be placed on the trajectory by adjusting
%               the angle between its segments (i.e., the final portion of
%               trajectory is not covered, in lack of information where the 
%               third point would be placed).
%
% See also INTERPTOUNIFORMSEGMENTS_FIXEDN

% Default parameters
if nargin < 5 
    symThresh = 0.05;
end
if nargin < 6 
    symTime = 1;
end
if nargin < 7
    decRound = 2;
end

% Make sure inputs are column vectors
sx = size(x); sy = size(y);
if sum(sx>1) > 1 || sum(sy>1) > 1
    error('x and y must be vectors, not matrices.')   
end
if sx(1)==1; x = x'; end
if sy(1)==1; y = y'; end

xy = [x,y];
localAngles = [];
trSegLens = sqrt(sum(diff(xy).^2,2)); % lengths of individual segments of input trajectory

length_to_shift = 0;    % distance by which the next point needs to be moved along trajectory in order
                        % to be stepSize arc length from the last points, starting
                        % from first point of segment that the next point will lie in.

% Function to compute angle in radians between vectors u and v
vAng =  @(u,v) acos(dot(u,v) /( norm(u)*norm(v)));

crs = 0; % current segment of trajectory, iterate over these
while crs < size(xy,1)-1
    
    % find segment of raw data in which the next point will lie,
    % based on cumulative arc length
    crs_old = crs;
    crs = crs + find(cumsum(trSegLens(crs+1:end)) > length_to_shift, 1);
    if isempty(crs)
        break;
    end
    
    % Reduce length_to_shift by total length of skipped segments
    length_to_shift = length_to_shift - sum(trSegLens(crs_old+1:crs-1));
    
    % points defining current segment
    p1 = xy(crs,:);
    p2 = xy(crs+1,:);
    
    % vector with same direction as current segment, unit length
    u = (p2-p1)/norm(p2-p1);
    
    base_point = p1;
    while 1 % switches to next raw segment when exceeding current segment's length (but remaining length of current seg is carried over)
        
        % where to put first vertex of two-segment fragment
        base_point = base_point + u * length_to_shift;
        
        % current segment length NOT exceeded
        if norm(base_point-p1) < norm(p2-p1)
            
            xy_temp = xy(crs:end, :);
            xy_temp(1,:) = base_point;
            % Fit two-segment fragment to trajctory
            % (i.e., interpolate two segments with chosen segment length)
            [interpCurve,~] = interpToUniformSegments_fixedn(xy_temp(:,1), xy_temp(:,2), segLen, 0 , symThresh, symTime, decRound, 3);
            % Compute angle between these segments
            if size(interpCurve,1) == 3
                localAngles(end+1) = vAng(interpCurve(2,:)-interpCurve(1,:), interpCurve(3,:)-interpCurve(2,:));
            end
            last_used_base_point = base_point;
            length_to_shift = stepSize;
            
        % current segment exceeded
        else
            
            if length_to_shift == stepSize
                length_to_shift = length_to_shift - norm(p2 - last_used_base_point);
            elseif length_to_shift < stepSize
                length_to_shift = length_to_shift - norm(p2 - p1);
            end
            break;
            
        end
        
    end
    
end

maxAngle = max(localAngles);

end