function [sys,x0,str,ts] = dynamic_snn1(t,x,u,flag)

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
sizes.NumContStates  = 14;
sizes.NumDiscStates  = 0;
sizes.NumOutputs     = 2;
sizes.NumInputs      = 2;
sizes.DirFeedthrough = 1;
sizes.NumSampleTimes = 1;

sys = simsizes(sizes);
x0  = [-0.3;0;0;0;0;0;0;0;0;0;0;0;0;0] ;
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

ep=0;
f_hat=x(13);
f_hatdots = x(14);
theta_hat=x(1:12);
p=zeros;
sig=zeros;
M=4;
r=sin(t);
bm=1;
%r=sin(t)-0.3*cos(7*t)+0.3*sint(11*t);
f=u(1);
u_in=u(2);
am=1;
g=1;


for i=1:M
p(i)=theta_hat((i+M))*f_hat+theta_hat(i+(2*M));
sig(i)=1/(1+exp(-p(i)));
theta_sig(i)=(theta_hat(i)*sig(i));
end
f_dot =f + u_in;
f_hatdot=(-am*f_hatdots)+am*f_hat+sum(theta_sig(:))+u_in;
for i=1:M
theta_hatdot(i)=g*(f_hat-f_hatdots)*sig(i);
theta_hatdot(i+M)=g*(f_hat-f_hatdots)*(theta_hat(i))*f_hat*(sig(i))^2*exp(-p(i));
theta_hatdot(i+(2*M))=g*(f_hat-f_hatdots)*(theta_hat(i))*(sig(i))^2*exp(-p(i));

end
ep=(f_hat-f_hatdots);
sys = [theta_hatdot';f_dot;f_hatdot]

% end mdlDerivatives
%
%=============================================================================
% mdlOutputs
% Return the block outputs.
%=============================================================================
%

function sys=mdlOutputs(t,x,u)

ep=0;
f_hat=x(13);
f_hatdots = x(14);
theta_hat=x(1:12);
p=zeros;
sig=zeros;
M=4;
r=sin(t);
bm=1;
%r=sin(t)-0.3*cos(7*t)+0.3*sint(11*t);
f=u(1);
u_in=u(2);
am=1;
g=1;


for i=1:M
p(i)=theta_hat((i+M))*f_hat+theta_hat(i+(2*M));
sig(i)=1/(1+exp(-p(i)));
theta_sig(i)=(theta_hat(i)*sig(i));
end
f_dot =f + u_in;
f_hatdot=(-am*f_hatdots)+am*f_hat+sum(theta_sig(:))+u_in;
for i=1:M
theta_hatdot(i)=g*(f_hat-f_hatdots)*sig(i);
theta_hatdot(i+M)=g*(f_hat-f_hatdots)*(theta_hat(i))*f_hat*(sig(i))^2*exp(-p(i));
theta_hatdot(i+(2*M))=g*(f_hat-f_hatdots)*(theta_hat(i))*(sig(i))^2*exp(-p(i));

end
ep=(f_hat-f_hatdots);
ep=(f_hat-f_hatdots);
sys=[sum(theta_sig(:));f_hatdots];
