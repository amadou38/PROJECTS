function pelt = parcours(elt, S)
  
  single1 = left(elt, S);
  single2 = top(elt, S);
  single3 = right(elt, S);
  single4 = bottom(elt, S);
  %pelt = cell(1, numel(single1)+numel(single2)+numel(single3)+numel(single4));
  pelt = {single1{1}, single2{1}, single3{1}, single4{1}, single1{2:end}, single2{2:end}, single3{2:end}, single4{2:end}};
end
