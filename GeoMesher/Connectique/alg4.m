function nedg = alg4(e2neigh,nelt)
% AUTEUR : Diallo Amadou, 28/09/2020
flag = false(nelt,1); nedg = 0;

for e = 1:nelt
    neloc = length(e2neigh{e});
    for eloc = 1:neloc
        tmp = e2neigh{e}(eloc,:);
        if (tmp ~= 0)
            if (flag(tmp) == ~true)      
                nedg = nedg + 1;
            end
        else
            nedg = nedg + 1;
        end
    end
    flag(e) = true;                 
end

end