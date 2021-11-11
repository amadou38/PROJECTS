%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mesh generator version 1.0                                            %
%                                                                       %
% Author: Amadou Diallo                                                 %
% Date: January 2021                                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all;
addpath('Connectique/');
addpath('mesh_files/');
addpath('../Project Laplace/benchmark/');
addpath('../Project Elasticity/benchmark/');
name = 'mesh1';     % name of the mesh
num_ind = 7;   % Corresponding face indices in the 'ind' cell array. NB: good from 1 to 13 
N = 16;        % N^2 corresponds to the total of principal elements on the mesh       
nbT = 2;
L = [1 1];     % DOM = [0 L(1)]x[0 L(2)]
ind = indices();    % function containing the default indices of principal elements
N = N*(max(ind{num_ind}(1,:))-1); % size of the box that contain elements 
[elements, Nod2e, vertices, e2edg, e2neigh, connect, freeedg, boundary,freenod] = generate(ind{num_ind},N,L); % function that generate the mesh

save_mat(name, elements, Nod2e, vertices, e2edg, e2neigh, connect, freeedg, boundary,freenod); % function that save the mesh
mesh = load('mesh1.mat');

%U = reshape(1:numel(elements), numel(elements), 1);
U = 10*rand(numel(elements),1);
figure,
plot_solution1(mesh,U);
filename = "solNum2x545x1";
readfile
tit = 'Approximate Solution 1';
figure,
plot_solution(mesh,tit,U);
filename = "solExa2x545x1";
readfile
tit = 'Exact Solution 1';
figure,
plot_solution(mesh,tit,U);
filename = "mesh1";
%convert
%[U,Uex] = VEM(mesh);
%[nelt,e2nod,nedg,e2dofs,e2edg,freeedg,freenod] = connexions(mesh);