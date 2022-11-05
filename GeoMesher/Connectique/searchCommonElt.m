function e0 = searchCommonElt(Nod2e,n1,n2,e)
% AUTEUR : Diallo Amadou, 28/09/2020

i0 = 0; j = 1;

while ((i0 == 0) && (j <= length(Nod2e{n1})))
    
    e1 = Nod2e{n1}(j);
    if (e1 ~= e)
        i = 0; k = 1;
        while ((i == 0) && (k <= length(Nod2e{n2})))
            e2 = Nod2e{n2}(k);
            if (e1 == e2)
                i = 1;
            else
                k = k + 1;
            end
        end
        if (i == 1)
            i0 = 1;
        else
            j = j + 1;
        end
    else
        j = j + 1;
    end
end
if (i0 == 1)
    if ((e1 == e2) && (e1 ~= e))
        e0 = e1;
    else
        e0 = 0;
    end
else
    e0 = 0;
end

end