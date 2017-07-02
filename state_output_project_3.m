function [sys,x0,str,ts] = state_output_project_3(t,x,u,flag)

switch flag,

  %%%%%%%%%%%%%%%%%%
  % Initialization %
  %%%%%%%%%%%%%%%%%%
  case 0,
    [sys,x0,str,ts]=mdlInitializeSizes;

  %%%%%%%%%%%%%%%
  % Derivatives %
  %%%%%%%%%%%%%%%
  case 1,
    sys=mdlDerivatives(t,x,u);

  %%%%%%%%%%%
  % Outputs %
  %%%%%%%%%%%
  case 3,
    sys=mdlOutputs(t,x,u);

  %%%%%%%%%%%%%%%%%%%
  % Unhandled flags %
  %%%%%%%%%%%%%%%%%%%
  case { 2, 4, 9 },
    sys = [];

  %%%%%%%%%%%%%%%%%%%%
  % Unexpected flags %
  %%%%%%%%%%%%%%%%%%%%
  otherwise
    error(['Unhandled flag = ',num2str(flag)]);

end
% end csfunc

%
%=============================================================================
% mdlInitializeSizes
% Return the sizes, initial conditions, and sample times for the S-function.
%=============================================================================
%
function [sys,x0,str,ts]=mdlInitializeSizes

% Initial Conditions
sizes = simsizes;
sizes.NumContStates  = 8;
sizes.NumDiscStates  = 0;
sizes.NumOutputs     = 6;
sizes.NumInputs      = 6;
sizes.DirFeedthrough = 1;
sizes.NumSampleTimes = 1;

sys = simsizes(sizes);
x0  = [0.1 0 0 0 0 0 0 0]';
str = [];
ts  = [0 0];

% end mdlInitializeSizes
%
%=============================================================================
% mdlDerivatives
% Return the derivatives for the continuous states.
%=============================================================================
%
function sys=mdlDerivatives(t,x,u)

% Input Signal

qd = u(1:2);         % Desired position
qdp = u(3:4);        % Differentiation of position
qdpp = u(5:6);       % Double differentiation of position

% System Parameters
m1 = 1;
m2 = 1;
a1 = 1;
a2 = 1;
g = 9.8;
am = 1;

% Controller gains
kp = 100;
kd = 14;

% Tracking Error Signals
e = qd - [x(1) x(2)]';
ep = qdp - [x(3) x(4)]';

% Inertial Matrix Elements
M11 = (m1+m2)*a1^2+m2*a2^2 + 2*m2*a1*a2*cos(x(2));
M12 = m2*a2^2 + m2*a1*a2*cos(x(2));
M22 = m2*a2^2;

% Non-Linear function Elements In System Dynamics (Coriolis/Centripetal Matrix + Gravitational Vector)
N1 = -m2*a1*a2*(2*x(3)*x(4)+x(4)^2)*sin(x(2));
N1 = N1 + (m1+m2)*g*a1*cos(x(1)) + m2*g*a2*cos(x(1) + x(2));
N2 = m2*a1*a2*x(3)^2*sin(x(2)) + m2*g*a2*cos(x(1) + x(2));

% Controller Design (PD)
s1= qdpp(1) + kd*ep(1) + kp*e(1);
s2= qdpp(2) + kd*ep(2) + kp*e(2);
tau1= M11*s1 + M12*s2 + N1;
tau2= M12*s1 + M22*s2 + N2;

det= M11*M22 - M12*M12;
MI11= M22/det;
MI12= -M12/det;
MI22= M11/det;

% Uncertainties & Friction Terms
F1 = 0.45*x(3);%0;%
F2 = 0.45*x(4);%0;%
taud1 = 0.4*sin(30*t);%0;%
taud2 = 0.4*sin(30*t);%0;%
% F1 = 0;%
% F2 = 0;%
% taud1 = 0;%
% taud2 = 0;%

% State Equations
xdot(1)= x(3);
xdot(2)= x(4);
% Friction Terms And Disturbances Added
xdot(3)= MI11*(tau1 - N1 - F1 - taud1) + MI12*(tau2 - N2 - F2 - taud2);
xdot(4)= MI12*(tau1 - N1 - F1 - taud1) + MI22*(tau2 - N2 - F2 - taud2);
xdot(5)= -am*(x(5)-x(3)) + MI11*(tau1 - N1) + MI12*(tau2 - N2);
xdot(6)= -am*(x(6)-x(4)) + MI12*(tau1 - N1) + MI22*(tau2 - N2);
xdot(7)= -am*(x(5)-x(3)) + MI11*(-F1 - taud1) + MI12*(-F2 - taud2);
xdot(8)= -am*(x(6)-x(4)) + MI12*(-F1 - taud1) + MI22*(-F2 - taud2);
sys = xdot';

% end mdlDerivatives
%
%=============================================================================
% mdlOutputs
% Return the block outputs.
%=============================================================================
%
function sys=mdlOutputs(t,x,u)
qd = u(1:2);

% Tracking Error For Postion
e= qd - [x(1) x(2)]';

% Output Equations
sys = [x(3);x(7);x(4);x(8);[x(5)-x(3);x(6)-x(4)]];
% end mdlOutputs