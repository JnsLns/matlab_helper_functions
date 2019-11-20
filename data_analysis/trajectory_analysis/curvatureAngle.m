function k = curvatureAngle(vertices)
% returns curvature of an ordered list of vertices that define a plane curve
% (m-by-2 matrix where m are vertices and columns are x/y coordinates).
% The returned curvature measure for a given vertex is the angle between
% the two bordering line segments, normalized by the sum of the segment's
% half lengths. The function yields nan for first and last data point.

% v are vertex coordinates

% in variable names the three vertices involved in curvature computation 
% for each individual vertice are denoted by a, b and c, corresponding to
% vertices v_i-1, v_i, and v_i+1

% vectors btw previous vertice and current vertice (v_i - v_i-1)
ab = diff(vertices);
ab = ab(1:end-1,:);
% vectors btw current vertice and next vertice (v_i+1 - v_i)
bc = diff(vertices);
bc = bc(2:end,:);

% angle between ab and bc
lambda = acos(sum(ab.*bc,2) ./ sqrt(sum(ab.^2,2) .* sum(bc.^2,2)));

% norm of ab
norm_ab = sqrt(sum(ab.^2,2));
% norm of bc
norm_bc = sqrt(sum(bc.^2,2));

% curvature is lambda normliazed by sum of halved lengths of the two
% surrounding edges
k = lambda;
k = [nan;k;nan];