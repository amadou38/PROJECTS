function [inElt, outElt, uniq] = other_elements(ind, pelt, A, Elt, n, m)

inElt = cell(1,n*m*4);
outElt = cell(1,n*m*(numel(pelt)-4));
nbre1 = 0; l1 = 0;
k = 1;
uniq = zeros(1,4*n*m);
p1 = 0;
y = 1:max(ind(2,:));
for j = 1:m
  x = 1:max(ind(1,:));
  for i = 1:n
    S = A(x,y);
    elt = Elt(k,:);
    pelt = parcours(elt, S);
    l2 = l1 + numel(pelt)-4;
    for l = 1:numel(pelt)-4
      outElt{l1+l} = pelt{4+l};
    end
    l1 = l2;
    p2 = 4*k;
    uniq(p1+1:p2) = [S(1,1) S(end,1) S(end,end) S(1,end)];
    p1 = p2;
    x = x + max(ind(1,:)) - 1;
    k = k + 1;
    nbre2 = nbre1 + 4;
    for count = 1:4
      inElt{nbre1+count} = pelt{count};
    end
    nbre1 = nbre2;
  end
  y = y + max(ind(2,:)) - 1; 
end

end