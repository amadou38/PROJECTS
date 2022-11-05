function Elt = princ_elements(ind, A, n, m)

% All principal elements
ni = max(ind(1,:)) - 1; %mj = max(ind(2,:));
nind = size(ind,2);
Elt = zeros(n*m,nind);
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

end