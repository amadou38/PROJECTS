function [U,Uex] = VEM(mesh)

% Main function
% AUTEUR : Diallo Amadou, 28/09/2020

addpath('../CodeAmadou/VEM_Laplace/');
%mesh.vertices(:,2) = 1.1*mesh.vertices(:,2);

[U, Uex] = vem_Laplace2(mesh);

end