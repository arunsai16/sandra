function [sys,x0,str,ts] = dyn_snn(t,x,u,flag)

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
sizes.NumOutputs     = 17;
sizes.NumInputs      = 1;
sizes.DirFeedthrough = 1;
sizes.NumSampleTimes = 1;

sys = simsizes(sizes);
x0  = [-0.3;0.1;0;0;0;0;0;0;0;0;0;0;0;0] ;
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

f_hat=x(13);
f_hatdots = x(14);
theta_hat=x(1:12);
theta_hatdot=zeros;
p=zeros;
sig=zeros;
xx=u;
M=4;
am=1;
g=[0.1;0.2;0.3;0.4;0.5;0.5;0.4;0.3;0.2;0.1;0.05;0.02;0.015;0.01];
f=((37.5*x(1)^2)+(50*x(1)^2*x(2)^2)+(3.5*x(2)));

for i=1:1:4
    
p(i)=theta_hat((i+M))*x(i)+theta_hat(i+(2*M));

end

for i=1:1:4

sig(i)=1/(1+exp(-p(i)));
theta_sig=(theta_hat(i)*sig(i));

end

f_dot = f + xx;
f_hatdot=(-am*f_hatdots)+(am*f_hat)+sum(theta_sig(:))+xx;
ep=(f_hat-f_hatdots);

for i=1:1:4
    
theta_hatdot(i)=g(i)*(ep)*sig(i);
theta_hatdot(i+M)=g(i+M)*(ep)*(theta_hat(i+M))*xx*(sig(i))^2*exp(-p(i));
theta_hatdot(i+(2*M))=g(i+(2*M))*(ep)*(theta_hat(i+(2*M)))*(sig(i))^2*exp(-p(i));

end

sys = [theta_hatdot';f_dot;f_hatdot]

% end mdlDerivatives
%
%=============================================================================
% mdlOutputs
% Return the block outputs.
%=============================================================================
%

function sys=mdlOutputs(t,x,u)

f_hat=x(13);
f_hatdots = x(14);
theta_hat=x(1:12);
p=zeros;
sig=zeros;
M=4;
xx=u;
f=((37.5*x(1)^2)+(50*x(1)^2*x(2)^2)+(3.5*x(2)));

for i=1:1:4

p(i)=x(i)*theta_hat(i+M)+theta_hat(i+(2*M));

end

for i=1:1:4

sig(i)=1/(1+exp(-p(i)));

end

ep=(f_hat-f_hatdots);

sys=[x;ep;f;f_hat];
