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
sizes.NumInputs      = 9;
sizes.DirFeedthrough = 1;
sizes.NumSampleTimes = 1;

sys = simsizes(sizes);
x0  = [0 0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.2 0.3 0.4]' ;
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

g=1;
am=1;
%inputs
tau1=u(1);
tau2=u(2);
N1=u(3);
N2=u(4);
MI11=u(5);
MI12=u(6);
qd_1=u(7);
xdot_3=u(9);
f_approx=x(1);
f_approx_state = x(2);
theta_hat=x(3:14);
M=4;
for i=1:M
    
p(i)=theta_hat((i+M))*qd_1+theta_hat(i+2*M);
sig(i)=1/(1+exp(-p(i)));%sigmoidal function with dynamic estimation
theta_sig(i)=theta_hat(i)*sig(i);
end

theta_sig=sum(theta_sig(:));
%as xdot_3 has the fault function it is estimated using 3 neurons with
%eight weights.
fault_sys=xdot_3;
f_esti_dot=-am*(f_approx_state-f_approx)+MI11*(tau1 - N1) + MI12*(tau2 - N2)+theta_sig;
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

qd_1=u(7);
q1=u(8);
xdot_3=u(9);
f_approx=x(1);
f_approx_state = x(2);
theta_hat=x(3:14);
M=4;

for i=1:M
    
p(i)=theta_hat(i+M)*qd_1+theta_hat(i+2*M);
sig(i)=1/(1+exp(-p(i)));
theta_sig(i)=theta_hat(i)*sig(i);
end

theta_sig=sum(theta_sig(:));

fault_system=xdot_3;
%f_esti_dot=-am*(f_esti_state-f_esti)+MI11*(tau1 - N1) + MI12*(tau2 - N2)+theta_sig;


ep=(f_approx-f_approx_state);

sys=[qd_1;q1;fault_system;theta_sig;ep];
