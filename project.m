function [sys,x0,str,ts] = project(t,x,u,flag)

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

sizes = simsizes;
sizes.NumContStates  = 8;
sizes.NumDiscStates  = 0;
sizes.NumOutputs     = 19;
sizes.NumInputs      = 0;
sizes.DirFeedthrough = 1;
sizes.NumSampleTimes = 1;

sys = simsizes(sizes);
x0  = [0 0 0.1 0 0 0 0 0]' ;
str = [];
ts  = [0 0];

% end mdlInitializeSizes
%
%=============================================================================
% mdlDerivatives
% Return the derivatives for the continuous states.
%=============================================================================
%
function sys=mdlDerivatives(t,x,~)
period = 6;
amp1 = 0.2;
amp2 = 0.2;
fact = 2*pi/period;
sinf = sin(fact*t);
cosf = cos(fact*t);
qd = [amp1*sinf amp2*cosf]';
qdp = fact*[amp1*cosf -amp2*sinf]';
qdpp = -fact^2*qd;
%pd computed-Torque control input
m1=1; m2=1; a1=1; a2=1; g=9.8; %arm parameters
kp=100; kv=14; %controller parameters
beta_fault=0;%fault profile
am=1;
%tracking errors
e=qd-[x(1) x(2)]';
ep=qdp-[x(3) x(4)]';
%computed inertia M(q) and nonlinear terms N(q,qdot)
M11=(m1+m2)*a1^2 +m2*a2^2 + 2*m2*a1*a2*cos(x(2));
M12=m2*a2^2 + m2*a1*a2*cos(x(2));
M22=m2*a2^2; 
N1=-m2*a1*a2*(2*x(3)*x(4)+ x(4)^2)*sin(x(2));
N1=N1+(m1+m2)*g*a1*cos(x(1))+m2*g*a2*cos(x(1)+x(2));
N2=m2*a1*a2*x(3)^2*sin(x(2))+m2*g*a2*cos(x(1)+x(2));
%pd CT control torques
s1=qdpp(1)+kv*ep(1)+kp*e(1);
s2=qdpp(2)+kv*ep(2)+kp*e(2);
tau1=M11*s1+M12*s2+N1;
tau2=M12*s1+M22*s2+N2;
%inversion of M(q) (for large values of n, use least-squares)
det=M11*M22-M12*M12;
MI11=M22/det;
MI12=-M12/det;
MI22=M11/det;
%uncertainity functions
Fd1=0.45*qd(1);
Fd2=0.45*qd(2);
Td_1=0.2*sin(30*t);
Td_2=0.2*sin(30*t);

%nominal solution
%xdot(3)=MI11*(-N1+tau1)+MI12*(-N2+tau2);
%xdot(4)=MI12*(-N1+tau1)+MI22*(-N2+tau2);
%equations when uncertainty is introduced
%given fault functions
fault_1 = 37.5*x(1)^2 + 50*x(1)^2*x(2)^2 + 3.5*x(2);
fault_2 = 0;
if t>10
    beta_fault=1;
end
%state equations
xdot(1)=x(3);
xdot(2)=x(4);
xdot(3)=MI11*(tau1-N1-Fd1-Td_1)+MI12*(tau2-N2-Fd2-Td_2)+beta_fault*fault_1;
xdot(4)=MI12*(tau1-N1-Fd1-Td_1)+MI22*(tau2-N2-Fd2-Td_2)+beta_fault*fault_2;
%equations to calculate threshold and residual
x3_hat=x(5);
x4_hat=x(6);
x3_hat_dot=-am*(x3_hat-x(3))+MI11*(tau1-N1)+MI12*(tau2-N2);
x4_hat_dot=-am*(x4_hat-x(4))+MI12*(tau1-N1)+MI22*(tau2-N2);
x_3_tilda_dot=-am*(-x3_hat+x(3))-MI11*(Fd1+Td_1)-MI12*(Fd2+Td_2);
x_4_tilda_dot=-am*(-x4_hat+x(4))-MI12*(Fd1+Td_1)-MI22*(Fd2+Td_2);

sys=[xdot';x3_hat_dot;x4_hat_dot;x_3_tilda_dot;x_4_tilda_dot];



% end mdlDerivatives
%
%=============================================================================
% mdlOutputs
% Return the block outputs.
%=============================================================================
%
function sys=mdlOutputs(t,x,~)
period = 6;
amp1 = 0.2;
amp2 = 0.2;
fact = 2*pi/period;
sinf = sin(fact*t);
cosf = cos(fact*t);
qd = [amp1*sinf amp2*cosf]';
qdp = fact*[amp1*cosf -amp2*sinf]';
qdpp = -fact^2*qd;
%pd computed-Torque control input
m1=1; m2=1; a1=1; a2=1; g=9.8; %arm parameters
kp=100; kv=20; %controller parameters
beta_fault=0;%fault profile
am=1;
%tracking errors
e=qd-[x(1) x(2)]';
ep=qdp-[x(3) x(4)]';
%computed inertia M9q) and nonlinear terms N(q,qdot)
M11=(m1+m2)*a1^2 +m2*a2^2 + 2*m2*a1*a2*cos(x(2));
M12=m2*a2^2 + m2*a1*a2*cos(x(2));
M22=m2*a2^2; 
N1=-m2*a1*a2*(2*x(3)*x(4)+ x(4)^2)*sin(x(2));
N1=N1+(m1+m2)*g*a1*cos(x(1))+m2*g*a2*cos(x(1)+x(2));
N2=m2*a1*a2*x(3)^2*sin(x(2))+m2*g*a2*cos(x(1)+x(2));
%pd CT control torques
s1=qdpp(1)+kv*ep(1)+kp*e(1);
s2=qdpp(2)+kv*ep(2)+kp*e(2);
tau1=M11*s1+M12*s2+N1;
tau2=M12*s1+M22*s2+N2;
%inversion of M(q) (for large values of n, use least-squares)
det=M11*M22-M12*M12;
MI11=M22/det;
MI12=-M12/det;
MI22=M11/det;
%equations to calculate threshold and residual
x_3_tilda=((1-exp(-am*t))/am)*(0.5*abs(MI11*qd(1)+MI12*qd(2))+0.3*abs(MI11+MI12));%threshold
x_4_tilda=((1-exp(-am*t))/am)*(0.5*abs(MI12*qd(1)+MI22*qd(2))+0.3*abs(MI12+MI22));%threshold
x_3_tilda_dot=x(7);
x_4_tilda_dot=x(8);
residue_1=abs(x_3_tilda_dot);%residual
residue_2=abs(x_4_tilda_dot);%residual
%uncertainity functions
Fd1=0.45*qd(1);
Fd2=0.45*qd(2);
Td_1=0.2*sin(30*t);
Td_2=0.2*sin(30*t);
%nominal solution
%xdot(3)=MI11*(-N1+tau1)+MI12*(-N2+tau2);
%xdot(4)=MI12*(-N1+tau1)+MI22*(-N2+tau2);
%given fault functions
fault_1 = 37.5*x(1)^2 + 50*x(1)^2*x(2)^2 + 3.5*x(2);
fault_2 = 0;
if t>10
    beta_fault=1;
end
%equations when uncertainty and are introduced
xdot(3)=MI11*(tau1-N1-Fd1-Td_1)+MI12*(tau2-N2-Fd2-Td_2)+beta_fault*fault_1;
xdot(4)=MI12*(tau1-N1-Fd1-Td_1)+MI22*(tau2-N2-Fd2-Td_2)+beta_fault*fault_2;
sys = [qd(1);x(1);qd(2);x(2);ep(1);ep(2);x_3_tilda;residue_1;x_4_tilda;residue_2;tau1;tau2;N1;N2;MI11;MI12;qd(1);x(1);xdot(3)];
% end mdlOutputs