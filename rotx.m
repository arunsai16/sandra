function [ T ] = rotx( alpha )
%Function provides homogeneous transformation matrix to rotate
%a co-ordinate system about x an angle of alpha degrees
T=[1 0 0 0;
   0 cosd(alpha) -sind(alpha) 0;
   0 sind(alpha) cosd(alpha) 0;
   0 0 0 1]; 
end
