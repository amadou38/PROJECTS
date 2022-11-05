function e2neigh = alg2(e2nod,Nod2e,nelt)
% AUTEUR : Diallo Amadou, 28/09/2020
e2neigh = cell(1,nelt);

for i = 1:nelt
    n1 = e2nod{i};
    neloc = length(n1);
    e2neigh{i} = zeros(neloc,1);
    for j = 1:neloc
        jp = j+1;
        if (j == neloc)
            jp = 1;
        end
        e2neigh{i}(j) = searchCommonElt(Nod2e,n1(j),n1(jp),i);
    end
end
end