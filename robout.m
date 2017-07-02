function[qd,e]=robout(t,x)
%compute desired trajectory
period=2;amp1=0.1; amp2=0.1;
fact=2*pi/period;
sinf=sin(fact*t);
cosf=cos(fact*t);
qd=[amp1*sinf amp2*cosf]';
%tracking errors
%e=qd-x(:,1:2);
t0=0;
tf=10;
x0=[0.1 0 0 0]';
[t,x]=ode23('robctl',t0,tf,x0);
[qd,e]=robout(t,x);
plot(t,e)