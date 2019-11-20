function dist = mousePointDist(window, point)

% [dist] = mousePointDist(window, point)
%
% Input arguments: window, vector[x,y] specifying the point (in pixels)
% to which the current pointer distance is to be calculated.
% Returns absolute distance in pixels.
%
% point can also be an n-by-2 row vector, giving one set of x,y coordinates
% in each row. In this case, the function returns a vector with n rows,
% each giving the current mouse distance from that point.


[x, y] = GetMouse(window);

dist = sqrt((x-point(:,1)).^2 + (y-point(:,2)).^2);
       

end