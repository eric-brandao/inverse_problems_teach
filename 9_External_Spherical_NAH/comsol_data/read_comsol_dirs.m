clc; clear all; close all;
%% loading and pre-processing pressure
%%%% Load real
filename = 'ls_150cm_30x60_real_1000.txt'; 
[nodes, triangles, pres_1k_r] = read_mesh_file(filename);
filename = 'ls_150cm_30x60_real_2000.txt'; 
[nodes, triangles, pres_2k_r] = read_mesh_file(filename);
filename = 'ls_150cm_30x60_real_3000.txt'; 
[nodes, triangles, pres_3k_r] = read_mesh_file(filename);

%%%% Load imag
filename = 'ls_150cm_30x60_imag_1000.txt'; 
[nodes, triangles, pres_1k_i] = read_mesh_file(filename);
filename = 'ls_150cm_30x60_imag_2000.txt'; 
[nodes, triangles, pres_2k_i] = read_mesh_file(filename);
filename = 'ls_150cm_30x60_imag_3000.txt'; 
[nodes, triangles, pres_3k_i] = read_mesh_file(filename);

pres = [pres_1k_r + 1j*[pres_1k_i; pres_1k_i(end)] , pres_2k_r + 1j*pres_2k_i, pres_3k_r + 1j*pres_3k_i];

%% Process nodes
[azimuth,elevation,r] = cart2sph(nodes(:,1),nodes(:,2),nodes(:,3));
r = 1.0*ones(size(r));
[nodes(:,1),nodes(:,2),nodes(:,3)] = sph2cart(azimuth,elevation,r);

scatter3(nodes(:,1),nodes(:,2),nodes(:,3), 10, pres_3k_i,'filled')
%%
save("sph_array_10x20.mat", "pres", "nodes", "triangles")