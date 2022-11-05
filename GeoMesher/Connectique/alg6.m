function connect = alg6(e2edg,nedg,nelt)
% AUTEUR : Diallo Amadou, 28/09/2020
connect = cell(1,nelt);
edgtagged = false(nedg,1);
for e = 1:nelt
    edgloc = e2edg{e};
    nedgloc = length(edgloc);
    for ee = 1:nedgloc
        if (edgtagged(edgloc(ee)) == false)
            connect{e}(ee) = edgloc(ee);
            edgtagged(edgloc(ee)) = true;
        else
            connect{e}(ee) = - edgloc(ee);
        end
    end
end

end