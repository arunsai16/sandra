function [sys,x0,str,ts] = static_rbf(t,x,u,flag)

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
sizes.NumContStates  = 12;
sizes.NumDiscStates  = 0;
sizes.NumOutputs     = 15;
sizes.NumInputs      = 1;
sizes.DirFeedthrough = 1;
sizes.NumSampleTimes = 1;

sys = simsizes(sizes);
x0  = [0.05;0.1;0.15;0.2;0.25;0.3;0.35;0.4;0.45;0.5;0.55;0.6] ;
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
func_hat=zeros;
for i=1:1:12
    
gamma=1;
theta_hat=x;
c=-1:(2/11):1;
sigmma=0.3;
xx=u;
numerator=(sin(2.5*xx)-(0.4*xx*(8+xx.^2)));
denominator=(0.5*(7+xx.^2));
func=numerator/denominator;
func_hat(i)=theta_hat(i)*exp(-(xx-c(i))^2/sigmma^2);

func_hat=sum(func_hat(:));
ep=(func-func_hat);
diff_func_hat=exp(-(xx-c).^2/sigmma.^2);

theta_hat_dot=gamma.*ep.*diff_func_hat;


end

sys = theta_hat_dot;

% end mdlDerivatives
%
%=============================================================================
% mdlOutputs
% Return the block outputs.
%=============================================================================
%

function sys=mdlOutputs(~,x,u)
ep=0;
func_hat=zeros;
for i=1:1:12
    

theta_hat=x;
c=-1:(2/11):1;
sigmma=0.3;
xx=u;
numerator=(sin(2.5*xx)-(0.4*xx*(8+xx.^2)));
denominator=(0.5*(7+xx.^2));
func=numerator/denominator;
func_hat(i)=theta_hat(i)*exp(-(xx-c(i)).^2/sigmma.^2);
func_hat=sum(func_hat(:));
%diff_func_hat=exp(-(xx-c).^2/sigmma.^2);
%theta_hat=gamma*(func-func_hat)*diff_func_hat;
%theta_hat_dot=gamma.*ep.*diff_func_hat;

ep=(func-func_hat);
end
sys=[x;ep;func;func_hat];

% end mdlOutputs





