function [sys,x0,str,ts] = state_output_project_2(t,x,u,flag)

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
sizes.NumOutputs     = 20;
sizes.NumInputs      = 6;
sizes.DirFeedthrough = 1;
sizes.NumSampleTimes = 1;

sys = simsizes(sizes);
x0  = [0 0 0.1 0 0 0 0 0]';
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

q1 = x(1);
q2 = x(2);
q1_dot = x(3);
q2_dot = x(4);

% System Parameters
m1 = 1;
m2 = 1;
a1 = 1;
a2 = 1;
g = 9.8;
theta1 = 0.45;
theta2 = 0.2;
am = 1;

% Controller gains
kp = 100;
kd = 14;

% Tracking Error Signals
e = qd - [q1 q2]';
ep = qdp - [q1_dot q2_dot]';

% Inertial Matrix Elements
M11 = (m1+m2)*a1^2+m2*a2^2 + 2*m2*a1*a2*cos(q2);
M12 = m2*a2^2 + m2*a1*a2*cos(q2);
M22 = m2*a2^2;

% Non-Linear function Elements In System Dynamics (Coriolis/Centripetal Matrix + Gravitational Vector)
N1 = -m2*a1*a2*(2*q1_dot*q2_dot+q2_dot^2)*sin(q2);
N1 = N1 + (m1+m2)*g*a1*cos(q1) + m2*g*a2*cos(q1 + q2);
N2 = m2*a1*a2*q1_dot^2*sin(q2) + m2*g*a2*cos(q1 + q2);

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
F1 = theta1*q1_dot;
F2 = theta1*q2_dot;
taud1 = theta2*sin(30*t);
taud2 = theta2*sin(30*t);
% F1 = 0;%
% F2 = 0;%
% taud1 = 0;%
% taud2 = 0;%

% Fault Terms
Fault1 = 37.5*q1^2 + 50*q1^2*q2^2 + 3.5*q2;
Fault2 = 0;
% Fault1 = 75*q1^2 + 100*q1^2*q2^2 + 7*q2 + 17;
% Fault2 = 50*q1*q2 + 12.5;
Beta = 0;

% Fault Profile Function
if t > 10
    Beta = 1;
end


% State Equations
xdot(1)= q1_dot;
xdot(2)= q2_dot;
% Friction Terms And Disturbances Added
xdot(3)= MI11*(tau1 - N1 - F1 - taud1) + MI12*(tau2 - N2 - F2 - taud2) + Beta*Fault1;
xdot(4)= MI12*(tau1 - N1 - F1 - taud1) + MI22*(tau2 - N2 - F2 - taud2) + Beta*Fault2;
% Estimation for the threshold design
xdot(5)= -am*(x(5) - x(3)) + MI11*(tau1 - N1) + MI12*(tau2 - N2);
xdot(6)= -am*(x(6) - x(4)) + MI12*(tau1 - N1) + MI22*(tau2 - N2);
% State estimation error; state equation
xdot(7)= (-am*(x(3) - x(5)) - MI11*(F1+taud1) - MI12*(F2+taud2));
xdot(8)= (-am*(x(4) - x(6)) - MI12*(F1+taud1) - MI22*(F2+taud2));
sys = xdot';

% end mdlDerivatives
%
%=============================================================================
% mdlOutputs
% Return the block outputs.
%=============================================================================
%
function sys=mdlOutputs(t,x,u)
qd = u(1:2);         % Desired position
qdp = u(3:4);        % Differentiation of position
qdpp = u(5:6);       % Double differentiation of position

q1 = x(1);
q2 = x(2);
q1_dot = x(3);
q2_dot = x(4);

% System Parameters
m1 = 1;
m2 = 1;
a1 = 1;
a2 = 1;
g = 9.8;
theta1 = 0.45;
theta2 = 0.2;
am = 1;

% Controller gains
kp = 100;
kd = 20;

% Tracking Error Signals
e = qd - [q1 q2]';
ep = qdp - [q1_dot q2_dot]';

% Inertial Matrix Elements
M11 = (m1+m2)*a1^2+m2*a2^2 + 2*m2*a1*a2*cos(q2);
M12 = m2*a2^2 + m2*a1*a2*cos(q2);
M22 = m2*a2^2;

% Non-Linear function Elements In System Dynamics (Coriolis/Centripetal Matrix + Gravitational Vector)
N1 = -m2*a1*a2*(2*q1_dot*q2_dot+q2_dot^2)*sin(q2);
N1 = N1 + (m1+m2)*g*a1*cos(q1) + m2*g*a2*cos(q1 + q2);
N2 = m2*a1*a2*q1_dot^2*sin(q2) + m2*g*a2*cos(q1 + q2);

% Controller Design (PD)
s1= qdpp(1) + kd*ep(1) + kp*e(1);
s2= qdpp(2) + kd*ep(2) + kp*e(2);
tau1= M11*s1 + M12*s2 + N1;
tau2= M12*s1 + M22*s2 + N2;

det= M11*M22 - M12*M12;
MI11= M22/det;
MI12= -M12/det;
MI22= M11/det;
% Tracking Error For Postion
e= qd - [q1 q2]';
% threshold1 = (((MI11)*(0.5*(q1_dot) + 0.3)) + ((MI12)*(0.5*(q2_dot) + 0.3)));
% threshold2 = (((MI12)*(0.5*(q1_dot) + 0.3)) + ((MI22)*(0.5*(q2_dot) + 0.3)));
% threshold_design1 = (threshold1)*((1-exp(-am*t))/am);
% threshold_design2 = (threshold2)*((1-exp(-am*t))/am);

threshold1 = 0.5*abs(MI11*q1_dot + MI12*q2_dot) + 0.3*abs(MI12 + MI11);
threshold2 = 0.5*abs(MI12*q1_dot + MI22*q2_dot) + 0.3*abs(MI22 + MI12);
threshold_design1 = (threshold1*(1-exp(-am*t))/am);
threshold_design2 = (threshold2*(1-exp(-am*t))/am);

residual1 = abs(x(7));
residual2 = abs(x(8));

F1 = theta1*q1_dot;
F2 = theta1*q2_dot;
taud1 = theta2*sin(30*t);
taud2 = theta2*sin(30*t);

% Fault Terms
Fault1 = 37.5*q1^2 + 50*q1^2*q2^2 + 3.5*q2;
Fault2 = 0;

% Fault1 = 75*q1^2 + 100*q1^2*q2^2 + 7*q2 + 17;
% Fault2 = 50*q1*q2 + 12.5;

Beta = 0;

% Fault Profile Function
if t > 10
    Beta = 1;
%     residual1 = abs(x(7));
%     residual2 = abs(x(8));
end

xdot(3)= MI11*(tau1 - N1 - F1 - taud1) + MI12*(tau2 - N2 - F2 - taud2) + Beta*Fault1
xdot(4)= MI12*(tau1 - N1 - F1 - taud1) + MI22*(tau2 - N2 - F2 - taud2) + Beta*Fault2

% Output Equations
sys = [qd;q1;q2;e;threshold_design1;residual1;threshold_design2;residual2;xdot(3);xdot(4);tau1;tau2;N1;N2;MI11;MI12;MI22;q1];
% end mdlOutputs