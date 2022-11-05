function [nelt,Nod2e,e2nod,nedg,e2dofs,e2edg,e2neigh,freeedg,freenod] = connexions(elements, vertices)
% AUTEUR : Diallo Amadou, 28/09/2020

nelt = length(elements);
nnod = size(vertices,1);
e2nod = elements;

Nod2e = alg1(e2nod,nnod,nelt);
e2neigh = alg2(e2nod,Nod2e,nelt);
nedg = alg4(e2neigh,nelt);
Verts = vertices;
[e2edg,nedg,freeedg,freenod] = alg5(e2neigh,e2nod,Verts,nelt,nnod,nedg);
e2dofs = alg6(e2edg,nedg,nelt);
e2dofs = e2dofs';
e2edg = e2edg';
e2neigh = e2neigh';
end