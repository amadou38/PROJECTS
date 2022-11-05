
clear; clc; close all;
addpath('Connectique/');
[X,Y] = meshgrid(0:.1:.8);
vertices = [X(:), Y(:)];
[n,m] = size(X);
A = reshape(1:n*m, n, m);
nx = n; ny = m;
Abord = A; Abord(2:end-1,2:end-1) = 0;
ind = [3 5 3 1 1; 1 3 5 4 2];                     % pentagon
%ind = [1 3 4 4 5 4 4 3 1 2; 1 2 1 2 3 4 5 4 5 3]; % stars
%ind = [2 4 5 5 4 2 1 1; 1 1 2 4 5 5 4 2];         % octogon
%ind = [2 3 4 5 5 4 2 1 1; 1 3 1 2 4 5 5 4 2];
%ind = [1 2 3 3 2 1; 1 1 2 3 3 2];                  % hexagon
%ind = [2 4 5 4 2 1 1; 1 2 3 4 5 4 2];               % heptagon
%ind = [1 2 3 4 5 5 4 5 5 4 3 2 1 1 2 1; 1 1 2 1 1 2 3 4 5 5 4 5 5 4 3 2];
%ind = [ 2 4 5 4 2 1; 1 1 3 5 5 3];
%ind = [3 5 5 3 1 1; 1 2 4 5 4 2];
%ind = [1 3 5 4 5 3 1 2; 1 2 1 3 5 4 5 3];
ni = max(ind(1,:)); mj = max(ind(2,:)); 
s = [1 1 ni ni; 1 mj mj 1];
ni = ni - 1; mj = mj - 1;
nind = size(ind,2);
n = (n-1)/ni; m = (m-1)/mj;
Elt = zeros(n*m,nind);
ni = max(ind(1,:)) - 1; 
nind = size(ind,2);
k2 = -ni;
for i = 1:n*m
    k1 = ni*mod(i-1,n);
    if k1 == 0
        k2 = k2 + ni;
    end
    for j = 1:nind
        Elt(i,j) = A(ind(1,j)+k1, ind(2,j)+k2);
    end
end

elts = cell(1,n*m);
for i = 1:n*m
    elts{i} = Elt(i,:);
end
size(vertices,1)
[A,elts,vertices,freenod] = organize(A,elts, vertices, nx, ny);
size(vertices,1)

for i = 1:numel(elts)
  for j = 1:numel(elts{i})
    Elt(i,j) = elts{i}(j);
  end
end
plus = n+1;
%plus = n+1+1;
%inElt = run(ind, Elt, A, n);
inElt = runright(ind, Elt, A, n);
elements = cell(1,n*m+plus);
elements = elts;
for i = 1:plus
  elements{n*m+i} = inElt{i};
end

name = 'mesh1';
save_mat(name, elements, vertices);
mesh = load('mesh_files/mesh1.mat');
%U = reshape(1:n*m, n*m, 1);
U = reshape(1:length(elements), length(elements), 1);
figure,
plot_solution(mesh,U);

