function [e2edg,nedg,freeedg,freenod] = alg5(e2neigh,e2nod,vtx,nelt,nnod,nedg)
% AUTEUR : Diallo Amadou, 28/09/2020
freeedg = false(nedg,1);
freenod = false(nnod,1);
nedg = 0;
e2edg = cell(1,nelt);
for e = 1:nelt
    nnodloc = length(e2nod{e});
    e2edg{e} = zeros(nnodloc,1);
end

for e = 1:nelt
    nodloc = e2nod{e};
    xt = vtx(nodloc,:);
    for eloc = 1:length(nodloc)
        if (e2edg{e}(eloc) ~= 0)
            continue;
        end
        if (e2neigh{e}(eloc) ~= 0)
            elocp = eloc + 1;
            if (eloc == length(nodloc))
                elocp = 1;
            end
            xmt = (xt(eloc,:) + xt(elocp,:))/2;
            nodlocneigh = e2nod{e2neigh{e}(eloc)};
            xl = vtx(nodlocneigh,:);
            for u =1:length(nodlocneigh) 
                up = u + 1;
                if (u == length(nodlocneigh))
                    up = 1;
                end
                xml = (xl(u,:) + xl(up,:))/2;
                if (norm(xmt-xml) <= 1e-18)
                    nedg = nedg + 1; 
                    e2edg{e}(eloc) = nedg;
                    e2edg{e2neigh{e}(eloc)}(u) = nedg;
                    break
                end
            end
        else
            nedg = nedg + 1;
            e2edg{e}(eloc) = nedg;
            freeedg(nedg) = true;
        end
    end    
end
% Dernier post-traitement
for e = 1:nelt
    neloc = length(e2nod{e});
    for eloc = 1:neloc
        if freeedg(e2edg{e}(eloc)) == true
            elocp = eloc + 1;
            if (eloc == neloc)
                elocp = 1;
            end
            freenod(e2nod{e}(eloc)) = true;
            freenod(e2nod{e}(elocp)) = true;
        end
    end
end

end