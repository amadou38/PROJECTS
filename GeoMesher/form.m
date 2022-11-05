function [inElt, outElt] = form(inElt,outElt,uniq,bounds,vertices)
  
  uniq = unique(uniq);
  uniq = uniq(~ismember(uniq,bounds));
  [~,Nod2e,~,~,~,~] = getIDs(inElt, vertices);
  noeud = cell(1,numel(uniq)+4);
  for i = 1:numel(uniq)
    intersect = Nod2e{uniq(i)};
    if numel(intersect)==2
      if uniq(i)>bounds(1) && uniq(i)<bounds(2)
        noeud{i}(1) = inElt{intersect(1)}(1);
        l = numel(inElt{intersect(2)}(2:end));
        noeud{i}(2:1+l) = inElt{intersect(2)}(2:end);
        noeud{i}(2+l:1+l+numel(inElt{intersect(1)}(3:end))) = inElt{intersect(1)}(3:end);
      elseif uniq(i)>bounds(2) && uniq(i)<bounds(3)
        noeud{i}(1) = inElt{intersect(1)}(1);
        l = numel(inElt{intersect(2)}(3:end));
        noeud{i}(2:1+l) = inElt{intersect(2)}(3:end);
        noeud{i}(2+l:1+l+numel(inElt{intersect(1)}(3:end))) = inElt{intersect(1)}(3:end);
      elseif uniq(i)>bounds(3) && uniq(i)<bounds(4)
        noeud{i}(1) = inElt{intersect(1)}(1);
        l = numel(inElt{intersect(2)}(3:end));
        noeud{i}(2:1+l) = inElt{intersect(2)}(3:end);
        noeud{i}(2+l) = inElt{intersect(2)}(1);
        noeud{i}(3+l:2+l+numel(inElt{intersect(1)}(3:end))) = inElt{intersect(1)}(3:end);
      else
        l1 = numel(inElt{intersect(1)}(3:end));
        noeud{i}(1:l1) = inElt{intersect(1)}(3:2+l1);
        l2 = numel(inElt{intersect(2)}(2:end));
        noeud{i}(1+l1:l1+l2) = inElt{intersect(2)}(2:end);
      end
    elseif numel(intersect)==4
      l1 = numel(inElt{intersect(1)}(3:end));
      noeud{i}(1:l1) = inElt{intersect(1)}(3:end);
      l1 = numel(inElt{intersect(1)}(3:end));
      noeud{i}(l1+1) = inElt{intersect(1)}(1);
      l2 = numel(inElt{intersect(2)}(1:end-1));
      noeud{i}(2+l1:l1+l2) = inElt{intersect(2)}(3:end);
      noeud{i}(1+l1+l2) = inElt{intersect(2)}(1);
      l2 = 1+l1+l2;
      l1 = numel(inElt{intersect(4)}(3:end));
      noeud{i}(l2+1:l1+l2) = inElt{intersect(4)}(3:2+l1);
      l2 = l1+l2;
      l1 = numel(inElt{intersect(3)}(3:end));
      noeud{i}(l2+1:l1+l2) = inElt{intersect(3)}(3:end);
    end
  end
  vect = zeros(1,numel(bounds));
  for i = 1:numel(bounds)
    vect(i) = Nod2e{bounds(i)};
  end
  for i = 1:numel(vect)
    noeud{numel(uniq)+i} = inElt{vect(i)};
  end
  inElt = cell(1,numel(noeud));
  for e = 1:numel(noeud)
      inElt{e} = unique(noeud{e},'stable');
  end
  
  n_outElt = 0; I = cell(1,2); %I{1}(1) = 0;
  for e = 1:numel(outElt)
      n_outElt = n_outElt + 1;
      v1 = [outElt{e}(1) outElt{e}(end)];
      for i = e+1:numel(outElt)
          v2 = [outElt{i}(1) outElt{i}(end)];
          if ismember(v1,v2) == true(1,2)
              n_outElt = n_outElt - 1;
              I{1}(end+1) = i;
              I{2}(end+1) = e;
          end
      end
  end
  rest_elts = I{1};
  noeud2 = cell(1,n_outElt);
  k = 0;
  for e = 1:numel(outElt)
      if ismember(rest_elts,e) == false
          k = k + 1;
          noeud2{k} = outElt{e};
      end
      v1 = [outElt{e}(1) outElt{e}(end)];
      for i = e+1:numel(outElt)
          v2 = [outElt{i}(1) outElt{i}(end)];
          if ismember(v1,v2) == true(1,2) 
            noeud2{k} = outElt{e};
            noeud2{k}(end+1:end+numel(outElt{i})-2) = outElt{i}(2:end-1);
          end
      end
  end
  if ~isempty(rest_elts)
    outElt = noeud2;
  end
end
