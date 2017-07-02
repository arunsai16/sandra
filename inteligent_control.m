function [sys,x0,str,ts] = inteligent_control(t,x,u,flag)
 
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
sizes.NumOutputs     = 1;
sizes.NumInputs      = 2;
sizes.DirFeedthrough = 1;
sizes.NumSampleTimes = 1;
 
sys = simsizes(sizes);
x0  = [ones(12,1);-0.3;1];
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
ip_x = u(1);
fx = u(2);
theta = x(1:12);
x_sysstate = x(13);
x_hat_sysstate = x(14);

gamma = 1;
am = 1;
 
 
    for i=1:4
        p(i) =  theta(i+a)*ip_x+theta(i+2*a);
    end
    for i = 1:4
        sigma(i)=1/(1+exp(-p(i)));
        fx_hat(i) = theta(i)*sigma(i);
    end
 
epsilon = x_sysstate-x_hat_sysstate;
 
x_dot = fx+ip_x;
x_hat_dot = -am*x_hat_sysstate+am*x_sysstate+sum(fx_hat(:))+ip_x;
    for i=1:a
         theta_dot(i) = gamma*epsilon*sigma(i);
         theta_dot(i+a)= gamma*epsilon*theta(i+a)*ip_x*sigma(i)^2*exp(-p(i));
         theta_dot(i+2*a) = gamma*epsilon*theta(i+2*a)*sigma(i)^2*exp(-p(i));
        
         
    end
 
 
sys =[theta_dot';x_dot;x_hat_dot];
 
% end mdlDerivatives
%
%=============================================================================
% mdlOutputs
% Return the block outputs.
%=============================================================================
%
function sys=mdlOutputs(t,x,u)
ip_x = u(1);
theta = x(1:12);

gamma = 1;
 
for i = 1:4
    p(i) =  theta(i+a)*ip_x+theta(i+2*a);
end
for i=1:4
sigma(i)=1/(1+exp(-p(i)));
fx_hat(i)=theta(i)*sigma(i);
end 
fx_hat = fx_hat';
 
sys = sum(fx_hat(:)) ;
 
% end mdlOutputs

