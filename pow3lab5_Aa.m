mu_0=4*pi*10^-7;
mu_rc1=150;
mu_rc2=2500;
A_c=80*10^-6;
N=10;
l_c=0.094;
p_c=0.01;
f=10*10^3:10:10*10^6;
R_hs1=(8.*f.*mu_rc1.*mu_0.*A_c.*(N.^2))./l_c;
R_es1=(pi.*(mu_rc1.*mu_0.*A_c.*f.*N).^2)./(2.*p_c.*l_c);
R_cs1=R_hs1+R_es1;
R_hs2=(8.*f.*mu_rc2.*mu_0.*A_c.*(N.^2))./l_c;
R_es2=(pi.*(mu_rc2.*mu_0.*A_c.*f.*N).^2)./(2.*p_c.*l_c);
R_cs2=R_hs2+R_es2;
figure(1)
loglog(f,R_hs1,'-c','LineWidth',2)
hold on
grid on
loglog(f,R_es1,'--m','LineWidth',2)
hold on
grid on
loglog(f,R_cs1,'--k','LineWidth',2)
hold on
grid on
xlabel('{\it f} (hz)','FontSize',12)
ylabel('{\it R} (?)','FontSize',12)
title('Analysis of equivalent series core resistance {\it R} vs. frequency {\it f}', 'FontSize', 12)
legend('R_h_s@?_r_c=150','R_e_s@?_r_c=150','R_c_s@?_r_c=150')
figure(2)
loglog(f,R_hs2,'-c','LineWidth',2)
hold on
grid on
loglog(f,R_es2,'--m','LineWidth',2)
hold on
grid on
loglog(f,R_cs2,'--k','LineWidth',2)
hold on
grid on
xlabel('{\it f} (hz)','FontSize',12)
ylabel('{\it R} (?)','FontSize',12)
title('Analysis of equivalent series core resistance {\it R} vs. frequency {\it f}', 'FontSize', 12)
legend('R_h_s@?_r_c=2500','R_e_s@?_r_c=2500','R_c_s@?_r_c=2500')