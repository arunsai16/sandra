function [ T ] = rotz( theta )
%Function provides homogeneous transformation matrix to rotate
%a co-ordinate system about z an angle of theta degrees

T=[cosd(theta) -sind(theta) 0 0;
   sind(theta) cosd(theta) 0 0;
   0 0 1 0;
   0 0 0 1];
end

