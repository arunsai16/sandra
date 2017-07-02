di=-5;
%intitial angle
df=40;
%final angle
tf=4;
%time taken for the joint to move from intial angle to final angle 
t=0*pi:pi/100:4;
a0=di;
a1=0;
a2=3*(df-di)/tf.^2;
a3=-2*(df-di)/tf.^3;
o = a0+a1*t+a2*t.^2+a3*t.^3;
subplot(3,1,1);
plot(t,o);
o1=a1+2*a2*t+(3*(a3*t.^2));
o2=2*a2+(6*(a3*t));
subplot(3,1,2);
plot(t,o1);
subplot(3,1,3);
plot(t,o2);

