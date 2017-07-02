function [sys,x0,str,ts] = dynamic_rbf(t,x,u,flag)

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
x0  = [0.05;0.1;0.15;0.2;0.25;0.3;0.35;0.4;0.45;0.5;0.55;0.6;0;0] ;
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

func_hat=x(13);
func_hat_dot_state = x(14);
for i=1:1:12
    
g=1;
theta=x(1:12);
am=1;
c=-1:(2/11):1;
sigmma=0.3;
xx=u;
numerator=(sin(2.5*xx)-(0.4*xx*(8+xx.^2)));
denominator=(0.5*(7+xx.^2));
func=numerator/denominator;
z=exp(-(xx-c).^2/sigmma.^2);
theta_z=(theta*z);
theta_z=sum(theta_z(:));
func_hat_dot=(-am*func_hat_dot_state)+(am*func_hat)+theta_z+xx;
ep=(func_hat-func_hat_dot_state);
theta_hat_dot=g.*ep.*z;
func_dot = func + xx;

end

sys = [theta_hat_dot';func_dot;func_hat_dot];

% end mdlDerivatives
%
%=============================================================================
% mdlOutputs
% Return the block outputs.
%=============================================================================
%

function sys=mdlOutputs(~,x,u)
ep=0;
func_hat=x(13);
func_hat_dot_state=x(14);
for i=1:1:12
    
g=1;
theta=x;
am=1;
c=-1:(2/11):1;
sigmma=0.3;
xx=u;
numerator=(sin(2.5*xx)-(0.4*xx*(8+xx.^2)));
denominator=(0.5*(7+xx.^2));
func=numerator/denominator;
z=exp(-(xx-c).^2/sigmma.^2);
theta_z=(theta*z);
theta_z=sum(theta_z(:));
func_hat_dot=(-am*func_hat_dot_state)+(am*func_hat)+theta_z+xx;
ep=(func_hat-func_hat_dot_state);
theta_hat_dot=g.*ep.*z;
func_dot=func+xx;
end
sys=[x;ep;func;func_hat];
% end mdlOutputs

