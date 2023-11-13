%% Código funcional e curto pra leitura de dados COMSOL - JAI 23

% Sistema de limpeza

clc; clear; close all; % limpar workspace e command window

%% Leitura arranjo 1
filename = 'ls_R0100cm_10x20_real.txt';
[nodes, triangles, p_real] = read_mesh_file_teste_joao(filename);

filename = 'ls_R0100cm_10x20_imag.txt';
[nodes, triangles, p_imag] = read_mesh_file_teste_joao(filename);
triangles = triangles -1;
pres = p_real + 1j*p_imag;

%%
scatter3(nodes(:,1),nodes(:,2),nodes(:,3), 10, 'filled')
%%
save("sph_array_10x20_100cm.mat", "pres", "nodes", "triangles")


%% Leitura arranjo 2
filename = 'ls_R025cm_10x20_real.txt';
[nodes, triangles, p_real] = read_mesh_file_teste_joao(filename);

filename = 'ls_R025cm_10x20_imag.txt';
[nodes, triangles, p_imag] = read_mesh_file_teste_joao(filename);
triangles = triangles -1;
pres = p_real + 1j*p_imag;

%%
scatter3(nodes(:,1),nodes(:,2),nodes(:,3), 10, 'filled')
%%
save("sph_array_10x20_25cm.mat", "pres", "nodes", "triangles")

%% Leitura arranjo 3
filename = 'ls_R0100cm_30x60_real.txt';
[nodes, triangles, p_real] = read_mesh_file_teste_joao(filename);

filename = 'ls_R0100cm_30x60_imag.txt';
[nodes, triangles, p_imag] = read_mesh_file_teste_joao(filename);
triangles = triangles -1;
pres = p_real + 1j*p_imag;

%%
scatter3(nodes(:,1),nodes(:,2),nodes(:,3), 10, 'filled')
%%
save("sph_array_30x60_100cm.mat", "pres", "nodes", "triangles")

%% Leitura arranjo 4
filename = 'ls_R025cm_30x60_real.txt';
[nodes, triangles, p_real] = read_mesh_file_teste_joao(filename);

filename = 'ls_R025cm_30x60_imag.txt';
[nodes, triangles, p_imag] = read_mesh_file_teste_joao(filename);
triangles = triangles -1;
pres = p_real + 1j*p_imag;

%%
scatter3(nodes(:,1),nodes(:,2),nodes(:,3), 10, 'filled')
%%
save("sph_array_30x60_25cm.mat", "pres", "nodes", "triangles")

%% Leitura arranjo 5
filename = 'ls_R0200cm_30x60_real.txt';
[nodes, triangles, p_real] = read_mesh_file_teste_joao(filename);

filename = 'ls_R0200cm_30x60_imag.txt';
[nodes, triangles, p_imag] = read_mesh_file_teste_joao(filename);
triangles = triangles -1;
pres = p_real + 1j*p_imag;

%%
scatter3(nodes(:,1),nodes(:,2),nodes(:,3), 10, 'filled')
%%
save("sph_array_30x60_200cm.mat", "pres", "nodes", "triangles")
