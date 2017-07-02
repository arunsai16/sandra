aj=1e-4
m=0.6
e=8.85e-12;
q=1.6e-19;
nd=1.12e10;
vd=0:0.1:2;
for vbi=[0.7 1 1.3]
    a=sqrt((e.*q.*nd)./(2.*vbi));
    cjo=aj*a;
    b=sqrt((1-(vd./vbi)).^m)
    cj=cjo./b
    plot(vd,cj)
end