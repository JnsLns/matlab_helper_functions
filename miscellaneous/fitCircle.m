function [x_m, y_m, r] = fitCircle(x,y)
% Fit circle through three 2d points
%
% x: three-element vector with x coordinates of three points
% y: three-element vector with y coordinates of three points
% x_m: x coordinate of fitted circle center
% y_m: y coordinate of fitted circle center
% r: radius of fitted circle

    % Get current vertices and vertices around it
    x1 = x(1);
    y1 = y(1);
    x2 = x(2);
    y2 = y(2);
    x3 = x(3);
    y3 = y(3);
            
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
    
end