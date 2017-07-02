function [sys,x0,str,ts] = project_3_fault_addition(t,x,u,flag)

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
sizes.NumContStates  = 11;
sizes.NumDiscStates  = 0;
sizes.NumOutputs     = 2;
sizes.NumInputs      = 16;
sizes.DirFeedthrough = 1;
sizes.NumSampleTimes = 1;

sys = simsizes(sizes);
x0  = [0 0 0 0 0 0 0 0 0 0 0]' ;
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
threshold_design1 = u(1);
residual1 = u(2);
threshold_design2 = u(3);
residual2 = u(4);
x_dot3 = u(5);
x_dot4 = u(6);
tau1 = u(7);
tau2 = u(8);
N1 = u(9);
N2 = u(10);
MI11 = u(11);
MI12 = u(12);
MI22 = u(13);
qd = u(14:15);
q1= u(16);
gamma = 1;
am = 1;

x_state = x(1);
x_hat_state = x(2);
theta_hat = x(3:11);
M = 3;

% D_function = 0;
% if residual1 > threshold1
%     D_function = q_tilde;
% end

for i = 1:M
    p(i) = theta_hat(i+M)*qd(1)+theta_hat(i+2*M);
    sigma(i) = 1/(1+exp(-p(i)));
    fx_hat(i) = theta_hat(i)*sigma(i);
end

fx_hat = sum(fx_hat(:));

x_sys_with_fault = x_dot3;
fault_estimation1 = -am*(x_hat_state - x_state) + MI11*(tau1 - N1) + MI12*(tau2 - N2) + fx_hat;

epsilon = (x_state - x_hat_state);

for i = 1:M
    theta_hat_dot(i) = gamma*epsilon*sigma(i);
    theta_hat_dot(i+M) = gamma*epsilon*theta_hat(i+M)*qd(1)*(sigma(i))^2*exp(-p(i));
    theta_hat_dot(i+2*M) = gamma*epsilon*theta_hat(i+2*M)*(sigma(i))^2*exp(-p(i));
end

sys = [x_sys_with_fault;fault_estimation1;theta_hat_dot'];

% end mdlDerivatives
%
%=============================================================================
% mdlOutputs
% Return the block outputs.
%=============================================================================
%
function sys=mdlOutputs(t,x,u)
threshold_design1 = u(1);
residual1 = u(2);
threshold_design2 = u(3);
residual2 = u(4);
x_dot3 = u(5);
x_dot4 = u(6);
tau1 = u(7);
tau2 = u(8);
N1 = u(9);
N2 = u(10);
MI11 = u(11);
MI12 = u(12);
MI22 = u(13);
qd = u(14:15);
q1= u(16);
gamma = 1;
am = 1;

x_state = x(1);
x_hat_state = x(2);
theta_hat = x(3:11);
M = 3;

% D_function = 0;
% if residual1 > threshold1
%     D_function = q_tilde;
% end

for i = 1:M
    p(i) = theta_hat(i+M)*qd(1)+theta_hat(i+2*M);
    sigma(i) = 1/(1+exp(-p(i)));
    fx_hat(i) = theta_hat(i)*sigma(i);
end

fx_hat = sum(fx_hat(:))

x_sys_with_fault = x_dot3;
fault_estimation1 = -am*(x_hat_state - x_state) + MI11*(tau1 - N1) + MI12*(tau2 - N2) + fx_hat;

epsilon = (x_state - x_hat_state);
sys =  [qd(1);q1];

% end mdlOutputs