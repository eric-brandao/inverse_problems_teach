clc; clear all; close all

load('sonah_pk.mat');
%% vetores k
index = 2; k = k0(index); 
%k_vec = k * dirs;
%% Matriz H

%%
Dx = 0.06; Dy = 0.06;
n_k = 60;
% (linspace)
kx = linspace(-pi/Dx, pi/Dx, n_k);
ky = linspace(-pi/Dy, pi/Dy, n_k);
% (k-space resolution) 
delta_kx = kx(2) - kx(1); delta_ky = ky(2) - ky(1);
% meshgrid
[kx_grid, ky_grid] = meshgrid(kx,ky);
% flatten
kx_f = kx_grid(:);
ky_f = ky_grid(:);

% form Kz
ke_norm = (kx_f.^2 + ky_f.^2).^0.5;
%kz_f = zeros(len(kx_f));
% propagating part
%idp = np.where(ke_norm <= k0)[0]
kz_f = sqrt(k.^2 - (kx_f.^2+ky_f.^2));
% evanescent part
%ide = np.where(ke_norm > k0)[0]
%kz_f[ide] = -1j*np.sqrt(kx_f[ide]**2+ky_f[ide]**2-k0**2)

% Kernel e Matriz
k_vec = [kx_f, ky_f, kz_f];
kappa = sqrt(delta_kx*delta_ky./(2*pi*k.^2));
fz_ref = 1.0 * sqrt(k./abs(kz_f));

% source plane
zp = -1.0*max([Dx,Dy]);
recs = [r_vecs(:,1), r_vecs(:,2), r_vecs(:,3)-zp];

% kernel e  matriz
H = transpose(repmat(fz_ref, 1, length(recs(:,1)))) .*...
    kappa .* exp(-1j * recs * transpose(k_vec));
%H = exp(-1j * recs * transpose(k_vec));
[U,s,V] = csvd(H);
picard(U,s,p(:,index));
%%
figure()
scatter3(kx_f, ky_f, real(kz_f))

%%
[lam_gcv,Gfun,rega_gcv] = gcv(U,s,p(:,index));
[lam_lc, rho ,eta,rega_lc] = l_curve(U,s,p(:,index));
[lam_ncp,dist,rega_ncp] = ncp(U,s,p(:,index));