function [ T ] = roty( phi )
%Function provides homogeneous transformation matrix to rotate
%a co-ordinate system about y an angle of phi degrees
T=[cosd(phi) 0 sind(phi) 0;
   0 1 0 0;
   -sind(phi) 0 cosd(phi) 0;
   0 0 0 1];

end
