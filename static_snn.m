function [sys,x0,str,ts] = static_snn(t,x,u,flag)

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


p=zeros;
sigmma=zeros;
theta_hat=zeros;
theta_hat_dot=zeros;
for i=1:1:4

M=4;
xx=u;
gamma=[0.1;0.2;0.3;0.4;0.5;0.6;0.7;0.8;0.9;1.0;1.1;1.2];
theta_hat(i)=x(i);
theta_hat(i+M)=x(i+M);
theta_hat(i+(2*M))=x(i+(2*M));
p(i)=(theta_hat(i+M))*xx+theta_hat(i+(2*M));
sigmma(i)=1/(1+exp(-p(i)));
numerator=(sin(2.5*xx)-(0.4*xx*(8+xx.^2)));
denominator=(0.5*(7+xx.^2));
func=numerator/denominator;
func_hat=theta_hat(i)*sigmma(i);
func_hat=sum(func_hat(:));
theta_hat_dot(i)=gamma(i)*(func-func_hat)*sigmma(i);
theta_hat_dot(i+M)=gamma(i+M)*(func-func_hat)*(theta_hat(i+M))*xx*(sigmma(i))^2*exp(-p(i));
theta_hat_dot(i+(2*M))=gamma(i+(2*M))*(func-func_hat)*(theta_hat(i+(2*M)))*(sigmma(i))^2*exp(-p(i));


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
theta_hat=zeros;
p=zeros;
sigmma=zeros;

for i=1:1:4

M=4;
xx=u;
theta_hat(i)=x(i);
theta_hat(i+M)=x(i+M);
theta_hat(i+(2*M))=x(i+(2*M));
p(i)=xx*theta_hat(i+M)+theta_hat(i+(2*M));
sigmma(i)=1/(1+exp(-p(i)));
numerator=(sin(2.5*xx)-(0.4*xx*(8+xx.^2)));
denominator=(0.5*(7+xx.^2));
func=numerator/denominator;
func_hat=theta_hat(i)*sigmma(i);
func_hat=sum(func_hat(:));
ep=(func-func_hat);

end

sys=[x;ep;func;func_hat];

% end mdlOutputs
