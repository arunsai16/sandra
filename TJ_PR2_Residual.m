
% IMPLEMENTATON OF FAULT FUNCTION AND RESIDUAL

function [sys,x0,str,ts] = TJ_PR2_Residual(t,x,u,flag)

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
sizes.NumOutputs     = 10;
sizes.NumInputs      = 0;
sizes.DirFeedthrough = 1;
sizes.NumSampleTimes = 1;

sys = simsizes(sizes);
x0  = [0 0.1 0 0 0 0 0 0]';
str = [];
ts  = [0 0];

% end mdlInitializeSizes
%
%=============================================================================
% mdlDerivatives
% Return the derivatives for the continuous states.
%=============================================================================
function sys=mdlDerivatives(t,x,u)

Time=6;
w=2*pi/Time;
qd=[0.2*sin(w*t) 0.2*cos(w*t)]';
qd_dot=w*[0.2*cos(w*t) -0.2*sin(w*t)]';
qd_ddot=-w^2*qd;

m1=0.5;m2=2;a1=1;a2=1;g=9.8;th_1=0.45;th_2=0.2;
kp=100;kd=20;

ep=qd-[x(1) x(2)]';
ep_dot=qd_dot-[x(3) x(4)]';

M11=(m1+m2)*a1^2+m2*a2^2+2*m2*a1*a2*cos(x(2));
M12=m2*a2^2+m2*a1*a2*cos(x(2));
M22=m2*a2^2;

N1=-m2*a1*a2*(2*x(3)*x(4)+x(4)^2)*sin(x(2));
N1=N1+(m1+m2)*g*a1*cos(x(1))+m2*g*a2*cos(x(1)+x(2));
N2=m2*a1*a2*x(3)^2*sin(x(2))+m2*g*a2*cos(x(1)+x(2));
N=[N1;N2];

s1=qd_ddot(1)+kd*ep_dot(1)+kp*ep(1);
s2=qd_ddot(2)+kd*ep_dot(2)+kp*ep(2);
T1=M11*s1+M12*s2+N1;
T2=M12*s1+M22*s2+N2;
T=[T1;T2];

det=M11*M22-M12*M12;
MI11=M22/det;
MI12=-M12/det;
MI22=M11/det;
MI=[MI11 MI12;MI12 MI22];

F=[(th_1*qd(1)) (th_1*qd(2))]';
Td=[(th_2*sin(30*t)) (th_2*sin(30*t))]';

U=MI*(-F-Td);

xdot(1)=x(3);
xdot(2)=x(4);
xdot(3)=MI11*(T1-N1-F(1)-Td(1))+MI12*(T2-N2-F(2)-Td(2));
xdot(4)=MI12*(T1-N1-F(1)-Td(1))+MI22*(T2-N2-F(2)-Td(2));
xdot=[xdot(1);xdot(2);xdot(3);xdot(4)];
am=1;

x3_hat=x(5);
x4_hat=x(6);
x3_tuda_dot=x(7);
x4_tuda_dot=x(8);
%ep_1=abs((x(3)-x3_hat));
%ep_2=abs((x(4)-x4_hat));
x3hat_dot=-am*(x3_hat-x(3))+MI11*(T1-N1)+MI12*(T2-N2);
x4hat_dot=-am*(x4_hat-x(4))+MI12*(T1-N1)+MI22*(T2-N2);
x3_tuda_dot=-am*(-x3_hat+x(3))-MI11*(F(1)+Td(1))-MI12*(F(2)+Td(2));
x4_tuda_dot=-am*(-x4_hat+x(4))-MI12*(F(1)+Td(1))-MI22*(F(2)+Td(2));

sys=[xdot;x3hat_dot;x4hat_dot;x3_tuda_dot;x4_tuda_dot];



% end mdlDerivatives
%
%=============================================================================
% mdlOutputs
% Return the block outputs.
%=============================================================================
%
function sys=mdlOutputs(t,x,u)

Time=6;
w=2*pi/Time;
qd=[0.2*sin(w*t) 0.2*cos(w*t)]';
qd_dot=w*[0.2*cos(w*t) -0.2*sin(w*t)]';
qd_ddot=-w^2*qd;


fault_fn=[((37.5*x(1)^2)+(50*x(1)^2*x(2)^2)+(3.5*x(2)));0];

ep=qd - [x(1) x(2)]';

am=1;

m1=0.5;m2=2;a1=1;a2=1;g=9.8;th_1=0.45;th_2=0.2;
kp=100;kd=20;

ep=qd-[x(1) x(2)]';
ep_dot=qd_dot-[x(3) x(4)]';

M11=(m1+m2)*a1^2+m2*a2^2+2*m2*a1*a2*cos(x(2));
M12=m2*a2^2+m2*a1*a2*cos(x(2));
M22=m2*a2^2;

N1=-m2*a1*a2*(2*x(3)*x(4)+x(4)^2)*sin(x(2));
N1=N1+(m1+m2)*g*a1*cos(x(1))+m2*g*a2*cos(x(1)+x(2));
N2=m2*a1*a2*x(3)^2*sin(x(2))+m2*g*a2*cos(x(1)+x(2));
N=[N1;N2];

s1=qd_ddot(1)+kd*ep_dot(1)+kp*ep(1);
s2=qd_ddot(2)+kd*ep_dot(2)+kp*ep(2);
T1=M11*s1+M12*s2+N1;
T2=M12*s1 + M22*s2 + N2;
T=[T1;T2];

det=M11*M22-M12*M12;
MI11=M22/det;
MI12=-M12/det;
MI22=M11/det;
MI=[MI11 MI12;MI12 MI22];

F=[(th_1*qd(1)) (th_1*qd(2))]';
Td=[(th_2*sin(30*t)) (th_2*sin(30*t))]';


x3tuda=((1-exp(-(am*t)))/am)*(0.5*abs(MI11*x(3)+MI12*x(4))+0.3*abs(MI11+MI12));
x4tuda=((1-exp(-(am*t)))/am)*(0.5*abs(MI12*x(3)+MI22*x(4))+0.3*abs(MI12+MI22));
x3_tuda_dot=x(7);
x4_tuda_dot=x(8);
residue1=abs(x3_tuda_dot);
residue2=abs(x4_tuda_dot);

sys=[qd(1);x(1);qd(2);x(2);ep(1);ep(2);x3tuda;residue1;x4tuda;residue2];
% end mdlOutputs

