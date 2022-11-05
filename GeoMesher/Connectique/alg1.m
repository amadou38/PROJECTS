function Nod2e = alg1(e2nod,nnod,nelt)
% AUTEUR : Diallo Amadou, 28/09/2020
Nod2e = cell(nnod,1);
Nbnod2e = zeros(nnod,1);

for i = 1 : nelt
    e = e2nod{i};
    for j = 1:length(e)
        Nod = e2nod{i}(j);
        nb = Nbnod2e(Nod) + 1;
        Nbnod2e(Nod) = nb;
        Nod2e{Nod}(nb) = i;
    end
end

end