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
sizes.NumOutputs     = 5;
sizes.NumInputs      = 1;
sizes.DirFeedthrough = 1;
sizes.NumSampleTimes = 1;

sys = simsizes(sizes);
x0  = [-0.3 0 0 0 0 0 0 0 0 0 0 0 0 0]' ;
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
%states
for xx=-1:1;
num=(sin(2.5*xx)-(0.4*xx*(8+xx.^2)));
den=(0.5*(7+xx.^2));
f=num/den;
end
am=1;
bm=2;
g=1;
r=u;
f_approx=x(1);
f_approx_state = x(2);
theta_hat=x(3:14);
M=4;
for i=1:M
    
p(i)=theta_hat((i+M))*x(1)+theta_hat(i+2*M);
sig(i)=1/(1+exp(-p(i)));%sigmoidal function with dynamic estimation
theta_sig(i)=theta_hat(i)*sig(i);
end

theta_sig=sum(theta_sig(:));
uu=-theta_sig-am*f_approx+bm*r;
f_sys=f
f_esti_dot=-am*(f_approx_state-f_approx)+theta_sig+uu;
ep=(f_approx-f_approx_state);

for i=1:M
   %neural network df/dt equations 
theta_hatdot(i)=g*(ep)*sig(i);
theta_hatdot(i+M)=g*(ep)*theta_hat(i+M)*qd_1*(sig(i))^2*exp(-p(i));
theta_hatdot(i+2*M)=g*(ep)*theta_hat(i+2*M)*(sig(i))^2*exp(-p(i));

end

sys = [fault_sys;f_esti_dot;theta_hatdot'];

% end mdlDerivatives
%
%=============================================================================
% mdlOutputs
% Return the block outputs.
%=============================================================================
%

function sys=mdlOutputs(~,x,u)

for xx=-1:1;
num=(sin(2.5*xx)-(0.4*xx*(8+xx.^2)));
den=(0.5*(7+xx.^2));
f=num/den;
end
am=1;
bm=2;
g=1;
r=u;
f_approx=x(1);
f_approx_state = x(2);
theta_hat=x(3:14);
M=4;
for i=1:M
    
p(i)=theta_hat((i+M))*x(1)+theta_hat(i+2*M);
sig(i)=1/(1+exp(-p(i)));%sigmoidal function with dynamic estimation
theta_sig(i)=theta_hat(i)*sig(i);
end

theta_sig=sum(theta_sig(:));
uu=-theta_sig-am*f_approx+bm*r;
f_sys=f;
f_esti_dot=-am*(f_approx_state-f_approx)+theta_sig+uu;
for i=1:M
   %neural network df/dt equations 
theta_hatdot(i)=g*(ep)*sig(i);
theta_hatdot(i+M)=g*(ep)*theta_hat(i+M)*qd_1*(sig(i))^2*exp(-p(i));
theta_hatdot(i+2*M)=g*(ep)*theta_hat(i+2*M)*(sig(i))^2*exp(-p(i));

end



ep=(f_approx-f_approx_state);

sys=[theta_sig;ep];
