function [outputCurve,arcSpanned] = interpToUniformSegments_fixedn(x,y,segLen,excessHandling,symThresh,symTimeThresh,roundToDec,fixedNumber)
% [outputCurve,arcSpanned] = interpToUniformSegments_fixedn(x,y,segLen,excessHandling,symThresh,symTimeThresh,roundToDec,fixedNumber)
%
% Interpolate 2-dimensional plane curve such that each segment in the
% resulting curve has the same length. The interpolation is done
% iteratively, progressing from first to last elements in x and y.
%
% outputCurve is an n-by-2 matrix, listing n x/y-coordinates, one set for
% for each interpolation vertex, in its rows.
%
% x and y are vectors holding x and y coordinates of the input curve's
% vertices.
%
% segLen is a scalar defining the desired segment length (same units as in
% arguments x,y).
%
% excessHandling (optional, default 0):
% if 0, any excess portion of the curve beyond the last interpolated
% segment that would result in a segment shorter than the desired length
% is discarded (ensuring that the output curve includes only segments of
% the desired length but shifting the coordinates of the last curve vertex);
% if 1, the final segment will be shorter than segLen and ends at the same
% point as the original curve;
% if 2, the direction of the final segment will be of length segLen and its
% direction will be defined by vector btw. most recently generated inter-
% polation vertex and the last point of the input curve (i.e., the final
% vertex' coordinates will differ from the input curve's end point but
% segment length will be uniform).
%
% symThresh (optional, default 0.15):
% This enables checking for eligible interpolation points using symbolic
% math for cases where an obtained numeric solution for such a point
% returned no eligiblity but was close. This is needed because with numeri-
% cally very close data points in the input curve, rounding errors can lead
% to impossibility of placing an interpolation point. This makes symbolic
% math checks necessary to avoid getting stuck in the computation loop.
% The meaning and scaling of the value symThresh is based on the core equa-
% tion used in this function namely
%               segLen^2 == sum((p2+l*(p3-p2)-p1).^2)
% which is solved (numerically or symbolically) for l to find a point
% p_interp that is both a euclidean distance of segLen from p1 (two-element
% vector) and lies on the curve segment between p2 and p3
%               p_interp = p2+(p3-p2)*l
% where l is the l obtained from the first equation. For a curve segment
% to be judged eligible for an interpolation point, l thus needs to be
% between 0 and 1. Rounding errors using the numerical method may lead to
% an eligible point being missed if analytically its l would be only
% slightly below one or only slighly above zero. Therefore, symbolic checks
% (guaranteed to find any eligible point) are done if it is found that l is
% up to symThresh above one or up to symThresh below zero, to make sure no
% placement was missed. Setting symThresh to 0 disables symbolic checks.
% Note that symbolic checks are computationally costly and heavily increase
% required processing time. smyThresh should therefore be lowered depending
% on how relevant precise solutions are (and based on whether or not
% computation gets stuck with the data at hand, which may happen when the
% last eligible point on the curve is missed).
%
% symTimeThresh (optional, default 1): time threshold in seconds for
% switching to symbolic math mode (slow) when, starting at the most recently
% generated interpolated point, no eligible interpolation point is found in
% in any of the tested segments within the time threshold.
% This may happen if symThresh is low and eligible points are missed due to
% rounding errors despite symThresh (may lead to infinite loops). Symbolic
% mode is then activated and search for eligible vertices on curve segments
% following the most recently generated point is done with exact
% computations. Once a point is found and stored, symbolic mode is
% deactivated again.
%
% roundToDec (optional, default -1):
% determines to how many decimals coordinate input is rounded prior to
% computing interpolation. If -1, no rounding takes place, otherwise
% roundToDec must be an integer >= 0 which defines the number of decimals
% to which x and y coordinates of input vertices are rounded. This
% determines up to which distance input vertices are considered to be equal;
% of multiple directly successive and equal (i.e., redundant) vertices,
% only the first one is taken into account when computing the interpolated
% trajectory, while the rest is removed beforehand. Rounding is not strictly
% needed since symThresh has been introduced, but initial rounding can
% lead to MUCH faster calculation (at the expense of some precision, of
% course).
%
% fixedNumber (optional, default inf):
% Causes function to return after a fixed number of vertices has been
% found.
%
%
% arcSpanned is an optional output vector that holds the arc length of the
% input curve that is spanned (i.e., lies between the two vertices of) each
% individual segment of the interpolated curve. If excessHandling==2 then
% the last value in that vector will be nan. 
%
% NOTE: The function contains code for incrementally plotting the
% interpolation process for debugging; can be enabled by setting doPlot=1
% at the outset of the function code.

% enable plot for debugging
doPlot = 0;


if nargin < 4
    excessHandling = 0;
end

if nargin < 5
    symThresh = 0.15;
end

if nargin < 6
    symTimeThresh = 1;
end

if nargin < 7
    roundToDec = -1;
end

if nargin < 8
    fixedNumber = inf;
end

% Round coordinates to desired number of decimals
if roundToDec >= 0
    x = round(x,roundToDec);
    y = round(y,roundToDec);
end

x = x(:); y = y(:); % make sure these are column vectors

% Remove no-move steps (redundant vertices) in curve
d = sqrt(sum(diff([x,y]).^2,2));
d = [1;d]; % never remove very first data point
x(d == 0,:) = [];
y(d == 0,:) = [];
x_backup = x;
y_backup = y;

% FOR DEBUGGING Compute output curve
if doPlot
    hFig = figure; plot(x,y); axis equal;
end

% First point of interpolated curve is same as that of input curve
outputCurve = [x(1) y(1)];
arcSpanned = [];

% FOR DEBUGGING
if doPlot
    hold on; h = plot(outputCurve(:,1),outputCurve(:,2),'marker','x');
end

% loop and generate one vertex in each iteration
done = 0;
exceededSymTimeThresh = 0;
foundVertices = 0;
while ~done && foundVertices < fixedNumber - 1       
    
    alwaysUseSymModeForCurVertex = 0;
    if exceededSymTimeThresh
        %disp('symmode'); % for debugging
        alwaysUseSymModeForCurVertex = 1;
        exceededSymTimeThresh = 0;
    end
    
    % FOR DEBUGGING
    if doPlot
        try
            hold on; set(h,'xdata',outputCurve(:,1),'ydata',outputCurve(:,2));
        end
        pause(0.1); drawnow;
    end
    
    % starting point: first vertex of input curve (the input data is
    % successively trimmed below, changing the coordinates assigned to p11
    % and p12 here; for later iterations p1 is the most recently generated
    % vertex).
    p11 = x(1);
    p12 = y(1);
    
    % loop over curve segments from current starting point, checking for each
    % whether it includes a point that is separated from current starting
    % point by segLen (i.e., an intersection with a circle around starting
    % point with radius segLen). If so, store coordinates of that point and
    % go to next new segment, if not, check next segment in input curve.
    tStartCurVert = tic;
    for curSegStart = 1:numel(x)-1
        
        if (toc(tStartCurVert) > symTimeThresh) && ~alwaysUseSymModeForCurVertex
            exceededSymTimeThresh = 1;
            break;
        end
        
        % FOR DEBUGGING
        if doPlot
            try
                delete(hCurSeg);
            end
        end
        
        % start point of segment currently under scrutiny
        p21 = x(curSegStart);
        p22 = y(curSegStart);
        % end point of segment currently under scrutiny
        p31 = x(curSegStart+1);
        p32 = y(curSegStart+1);
        
        % skip segments of zero length
        if all([p31 p32]==[p21 p22])
            continue
        end
            
        % FOR DEBUGGING
        if doPlot
            hCurSeg = plot([p21 p31],[p22 p32],'marker','o');
        end
        
        % For the vector line equation  p_interp = p2+(p3-p2)*l ,
        % find 0<l<=1 such that ||p_interp-p1||=segLen.
        % Numeric approach (non-negative solution for l):
        if ~alwaysUseSymModeForCurVertex
            l = ...
                (p11*p31 - p12*p22 - p11*p21 + p12*p32 - p21*p31 - p22*p32 + ...
                (segLen^2*p21^2 - 2*segLen^2*p21*p31 + segLen^2*p22^2 - 2*segLen^2*p22*p32 + ...
                segLen^2*p31^2 + segLen^2*p32^2 - p11^2*p22^2 + 2*p11^2*p22*p32 - ...
                p11^2*p32^2 + 2*p11*p12*p21*p22 - 2*p11*p12*p21*p32 - ...
                2*p11*p12*p22*p31 + 2*p11*p12*p31*p32 - 2*p11*p21*p22*p32 + ...
                2*p11*p21*p32^2 + 2*p11*p22^2*p31 - 2*p11*p22*p31*p32 - ...
                p12^2*p21^2 + 2*p12^2*p21*p31 - p12^2*p31^2 + 2*p12*p21^2*p32 -...
                2*p12*p21*p22*p31 - 2*p12*p21*p31*p32 + 2*p12*p22*p31^2 - ...
                p21^2*p32^2 + 2*p21*p22*p31*p32 - p22^2*p31^2)^(1/2) + ...
                p21^2 + p22^2)/(p21^2 - 2*p21*p31 + p22^2 - 2*p22*p32 + p31^2 + p32^2);
        end
        
        % if l is just above one or just below zero, check exact value
        % using symbolic math:
        symbolicMode = 0;
        if  alwaysUseSymModeForCurVertex || ((l>=1-symThresh && l<=1+symThresh) || (l>=0-symThresh && l<=0+symThresh)) || ~isreal(l)
            symbolicMode = 1;            
            syms l;
            p1 = [p11 p12];
            p2 = [p21 p22];
            p3 = [p31 p32];
            r = segLen;
            eqn1 = r^2 == sum((p2+l*(p3-p2)-p1).^2);
            l = vpasolve(eqn1,l);            
        end
        
        % FOR DEBUGGING
        if doPlot
            l_tmp = double(l(logical(l > 0)));
            if ~isempty(l_tmp )
                try
                    delete(hPNew)
                end
                pNew = [p21 p22] + ([p31 p32]-[p21 p22]) * l_tmp;
                %d = norm(([p21 p22] + ([p31 p32]-[p21 p22]) * l_tmp)-[p11 p12]);
                %disp(['distance =' num2str(double(d))]);
                %disp(['l =' num2str(double(l_tmp))]);
                hPNew = plot(pNew(1),pNew(2),'marker','+');
            end
        end
        
        % if eligible l has been found, compute and store p_interp and leave
        % this loop (starting new seach from the new point); if not, check
        % next curve segment; do separate checks for eligibility of l, de-
        % pending on symbolic/numeric mode.
        if  (~symbolicMode && l <= 1 && l > 0) || ...
                (symbolicMode && any(logical(l <= 1) & logical(l > 0)))
            
            % get positive solution for l and convert to double
            if symbolicMode
                l = double(l(logical(l <= 1) & logical(l > 0)));
            end
            
            % append newly generated vertex (p_interp) to new curve
            outputCurve(end+1,:) = [p21, p22] + ([p31, p32]-[p21, p22]) * l;            
            % get arc length of input curve spanned by the current segment
            arcSpanned(end+1,1) = sum(sqrt(diff([x(1:curSegStart);outputCurve(end,1)]).^2 + diff([y(1:curSegStart);outputCurve(end,2)]).^2));
            % Remove segments from input data that are spanned by new segment            
            % of interpolated curve and replace first point with last
            % generated new point (this will be used as starting point
            % for next search for p_interp)                                                            
            x(1:curSegStart-1) = [];
            y(1:curSegStart-1) = [];
            x(1) = outputCurve(end,1);
            y(1) = outputCurve(end,2);
            % go to search for next new vertex
            break;
            % in case no eligible l is found for current segment and we are
            % at the final segment of the input curve, either
            % (excessHandling = 0) discard remaindor of input trajectory
            % (i.e., interpolated curve will end earlier) or
            % (excessHandling = 1) place final p_interp at position of last
            % point of input data (i.e., final segment of interpolated curve
            % will be shorter than segLen), or (excessHandling = 2) add
            % final segment of length segLen whose direction is determined
            % by vector between last added p_interp and final point of
            % input data.
        elseif curSegStart == numel(x)-1
            if excessHandling == 1
                outputCurve(end+1,1) = x(end);
                outputCurve(end,2) = y(end);
                arcSpanned(end+1,1) = sum(sqrt(diff([x(1:curSegStart);outputCurve(end,1)]).^2 + diff([y(1:curSegStart);outputCurve(end,2)]).^2));
            elseif excessHandling == 2
                outputCurve(end+1,:) = [p11, p12] + (([x(end),y(end)]-[p11,p12])/norm([x(end),y(end)]-[p11,p12]))*segLen;
                arcSpanned(end+1,1) = nan;
            end
            done = 1;
            break;
        end
        
    end
  
    foundVertices = foundVertices + 1;
    
end

% FOR DEBUGGING
if doPlot
    try
        close(hFig)
    end
end

end