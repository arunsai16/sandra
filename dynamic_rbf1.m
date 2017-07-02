function [sys,x0,str,ts] = dynamic_rbf1(t,x,u,flag)

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
function sys=mdlDerivatives(~,x,u)
f_hat=x(13);
f_hatdots = x(14);
f=u(1);
u_in=u(2);
for i=1:1:12
    
g=1;
theta=x(1:12);
am=1;
c=-1:(2/11):1;
sig=0.3;

z=exp(-(f_hat-c).^2/sig.^2);
theta_z=(theta*z);
theta_z=sum(theta_z(:));
f_hatdot=(-am*f_hatdots)+(am*f_hat)+theta_z+u_in;
ep=(f_hat-f_hatdots);
theta_hatdot=g.*ep.*z;
f_dot = f + u_in;

end

sys = [theta_hatdot';f_dot;f_hatdot];

% end mdlDerivatives
%
%=============================================================================
% mdlOutputs
% Return the block outputs.
%=============================================================================
%

function sys=mdlOutputs(~,x,u)
f_hat=x(13);
f_hatdots = x(14);
f=u(1);
u_in=u(2);
for i=1:1:12
    
g=1;
theta=x(1:12);
am=1;
c=-1:(2/11):1;
sig=0.3;

z=exp(-(f_hat-c).^2/sig.^2);
theta_z=(theta*z);
theta_z=sum(theta_z(:));
f_hatdot=(-am*f_hatdots)+(am*f_hat)+theta_z+u_in;
ep=(f_hat-f_hatdots)
theta_hatdot=g.*ep.*z;
f_dot = f + u_in;


end
sys=[theta_z;f_hatdots];
% end mdlOutputs

