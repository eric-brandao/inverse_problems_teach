clc; clear all; close all

load('ang_spk2.mat');
%% vetores k
index = 2; k = k0(index); 
k_vec = k * dirs;
%% Matriz H
H = exp(-1j*r_vecs * transpose(k_vec));
[U,s,V] = csvd(H);
picard(U,s,p(:,index));
%%
[lam_gcv,Gfun,rega_gcv] = gcv(U,s,p(:,index));
[lam_lc, rho ,eta,rega_lc] = l_curve(U,s,p(:,index));
%[lam_ncp,dist,rega_ncp] = ncp(U,s,p(:,index));

%%
[x_lambda,rho,eta] = tikhonov(U,s,V,p(:,index),lam_lc);