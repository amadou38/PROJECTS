function single = top(elt,S)
  
  %for i = 1:numel(elt)
    single = cell(1,1);
    [I,V1] = find(ismember(elt,S(:,1)));
    if I ~= 0
      single{1}(1) = elt(V1(end));
    end
    if ~ismember(S(end,1),elt)
      single{1}(2) = S(end,1);
    end
    [I,V2] = find(ismember(elt,S(end,:)));
    if I ~= 0
      single{1}(3) = elt(V2(1));
    end
    for k = 1:numel(elt(V1(end)+1:V2(1)))-1
      single{1}(3+k) = elt(V2(1)-k);
    end
    
    if abs(V2(1)-V2(end)) > 1        % Verify if indices form an element
      single{2} = elt(V2(1):V2(end));
    end
  %end
  
end
