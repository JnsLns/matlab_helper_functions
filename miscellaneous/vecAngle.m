function ang = vecAngle(u,v)
% Angle in radians between vectors u and v
ang = acos(dot(u,v) /( norm(u)*norm(v)));
end