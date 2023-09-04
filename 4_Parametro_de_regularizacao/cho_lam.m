clc; clear all; close all

load('gravity_64.mat');

[U,s,V] = csvd(A);

[x_delta,lam_dp] = discrep(U,s,V,gxln', norm(n));
[lam_gcv,Gfun,rega_gcv] = gcv(U,s,gxln');
[lam_lc, rho ,eta,rega_lc] = l_curve(U,s,gxln');
[lam_ncp,dist,rega_ncp] = ncp(U,s,gxln');

%%
[x_dp,rho,eta] = tikhonov(U,s,V,gxln',lam_dp);
[x_gcv,rho,eta] = tikhonov(U,s,V,gxln',lam_gcv);
[x_lc,rho,eta] = tikhonov(U,s,V,gxln',lam_lc);
[x_ncp,rho,eta] = tikhonov(U,s,V,gxln',lam_ncp);

figure()
plot(x, rhox); hold on;
plot(x, x_dp); hold on;
plot(x, x_gcv); hold on;
plot(x, x_lc); hold on;
plot(x, x_ncp); hold on;
legend('Ref', 'DP', 'GCV', 'L-c', 'NCP')