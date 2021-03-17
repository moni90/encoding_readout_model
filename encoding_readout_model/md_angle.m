function theta = md_angle(v,w)

% v and w are n x 1 vectors
% theta is the angle beetween these vectors

% theta = acos((v'*w)/(sqrt(v'*v)*sqrt(w'*w)));

theta = acos(dot(v,w)/(norm(v)*norm(w)));
if theta > pi && theta <= 2*pi 
   theta = theta-pi;
end

if theta > pi/2
    theta = pi-theta;
end

% theta = atan2(w(1)*v(2)-w(2)*v(1),dot(v,w));
% if abs(theta) > pi/2
%     theta = -(pi-theta);
% end
