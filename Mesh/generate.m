function [elements, Nod2e, vertices, e2edg, e2neigh, connect, freeedg, boundary,freenod] = generate(ind,N,L)
% Function that generate 2D meshes from the indices of the given element by
% the user or by using the default elements that are defined  in indices.m
%
% SYNOPSIS : [elements, vertices, boundary] = generate(ind,N,L)
%
% INPUT : ind: matrix 2xn containing the coordinates of all nodes that form 
%              the first principal element
%           N: size of the box that contain elements
%           L: DOM = [0 L(1)]x[0 L(2)]
%
% OUTPUT : elements: cell array [cell(nelt,1)] containing nodes of each elements
%          vertices: coordinates (Nx2) of each nodes 
%          boundary: all nodes on the boundary
%
% EXAMPLE : 
%           ind = [4 5 7 5 4 3 1 3; 1 3 4 5 7 5 4 3];   
%            N = 10; 
%            N = N*(max(ind(1,:))-1);
%            L = [1.5 1] 
%
%
h1 = (L(1)-0)/N;
h2 = (L(2)-0)/N;
[X,Y] = meshgrid(0:h1:L(1),0:h2:L(2));
vertices = [X(:), Y(:)];
[n,m] = size(X);
A = reshape(1:n*m, n, m);
nx = n; ny = m;
ni = max(ind(1,:)); mj = max(ind(2,:)); 
ni = ni - 1; mj = mj - 1;
n = fix((n-1)/ni); 
m = fix((m-1)/mj);

% First principal element
x = 1:max(ind(1,:)); y = 1:max(ind(2,:));
S = A(x,y);
elt = zeros(1,numel(ind(1,:)));
for j = 1:numel(ind(1,:))
    elt(j) = S(ind(1,j), ind(2,j));
end
pelt = parcours(elt, S);

Matrix_Elt = princ_elements(ind, A, n, m);

[inElt, outElt, uniq] = other_elements(ind, pelt, A, Matrix_Elt, n, m);

% A revoir
bounds = [1 nx (nx-1)*ny+1 nx*ny];
[inElt, outElt] = form(inElt,outElt,uniq,bounds,vertices);

elts = cell(1,n*m);
for i = 1:n*m
    elts{i} = Matrix_Elt(i,end:-1:1);
    this = elts{i};
    elts{i}(1) = this(end);
    elts{i}(2:end) = this(1:end-1);
end

for i = 1:numel(inElt)
    this = inElt{i};
    inElt{i} = this(end:-1:1);
end
% for i = 1:numel(outElt)
%     this = outElt{i};
%     outElt{i} = this(end:-1:1);
% end
elements = {elts{:} inElt{:} outElt{:}};
%elements = elts;
fprintf('Number of vertices before the reorganization: %d\n',size(vertices,1));
[elements,vertices,boundary] = organize(elements, vertices);
fprintf('Number of vertices after the reorganization: %d\n',size(vertices,1));   
fprintf('Number of elements: %d\n',numel(elements));

elements = elements';
first_ind = zeros(numel(elements),1);
for i = 1:numel(elements)
    elements{i} = elements{i}';
    first_ind(i) = min(elements{i});
end
[X,index] = sort(first_ind);
elements1 = elements;
for i = 1:numel(X)
    elements{i} = elements1{index(i)};
end

[~,Nod2e,~,~,connect,e2edg,e2neigh,freeedg,freenod] = connexions(elements, vertices);
freeedg = unique(find(freeedg));
freenod = unique(find(~freenod));

end

