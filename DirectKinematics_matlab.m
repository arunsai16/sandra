%
%  This is a Matlab program to calculate direct kinematics
%

% clean up previous workspace
clear all;
close all;
clc;

% define symbolic variables
syms d1 L1 L2 Lt t1 t2 t3 t4 t5 r1 r2 r3 r4 r5 r6;

% define the degree of freedom. For our robot N=5. For Puma N=6
N=5;
N=N+1;

% 
% Define DH table
% Insert your DH table here. The dimension of the DH table is (N+1)x4
% In this table the last row (N+!) relates the tool
% frame (6 in our case) to the Nth frame (5th frame in our case).  
% Remeber the tool frame (N+!) has the same orientation as the Nth frame
% Therefore, alpha(i-1) and theta(i-1) is zero. Also the translation is in
% the z-direction (either positive or negative depending on your coordinate
% system. Therefore, the only entry will be d(i) in the last row.
%
dh = [180,0,-d1,t1;
      -90,0,0,t2;  
      0,L1,0,t3;
      0,L2,0,t4-pi/2;
      90,0,0,t5;
      0,0,-Lt,0;];

% Find T01 to T(i-1)i
for i=1:N
    % 
    % 1st row of T
    %
    a11 = cos(dh(i,4));  
    a12 = -sin(dh(i,4));  
    a13 = 0;  
    a14 = dh(i,2);
    %
    % 2nd row of T
    %
    a21 = sin(dh(i,4))*cos(dh(i,1));
    a22 = cos(dh(i,4))*cos(dh(i,1));
    a23 = -sin(dh(i,1));
    a24 = -sin(dh(i,1))*dh(i,3);
    %
    % 3rd row of T
    %
    a31 = sin(dh(i,4))*sin(dh(i,1));
    a32 = cos(dh(i,4))*sin(dh(i,1));
    a33 = cos(dh(i,1));
    a34 = cos(dh(i,1))*dh(i,3);
    %
    % 4th row of T
    %
    a41 = 0;
    a42 = 0;
    a43 = 0;
    a44 = 1;
    
    % 
    % T for current index
    %
    t(:,:,i) = [a11, a12, a13, a14;
                a21, a22, a23, a24;
                a31, a32, a33, a34;
                a41, a42, a43, a44;
               ];
end

% calculate t05
t05 = t(:,:,1)*t(:,:,2)*t(:,:,3)*t(:,:,4)*t(:,:,5);

% 
% calculate T06. The method used in this m-file is to multiply T05 and T56. 
% 
t06 = t05*t(:,:,6);

% Simplify t0 and t06
t05 = simplify(t06);
t06 = simplify(t06);

% we can then fill in symbolic values with numbers, if we so desire*
% just uncomment the code below to try this

 d1 = 300; 
 L1 = 250; 
 L2 = 160;
 Lt = 179;
 t1=30;
 t2=20;
 t3=10;
 t4=-50;
 t5=10;

subs(t05);
subs(t06);
t;