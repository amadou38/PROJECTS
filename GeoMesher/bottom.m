function single = bottom(elt,S)
  
  %for i = 1:numel(elt)
    single = cell(1,1);
    [I,V1] = find(ismember(elt,S(:,end)));
    if I ~= 0
      single{1}(1) = elt(V1(end));
    end
    if ~ismember(S(1,end),elt)
      single{1}(2) = S(1,end);
    end
    [I,V1] = find(ismember(elt,S(1,:)));
    if I ~= 0
      single{1}(3) = elt(V1(1));
    end
    for l = 1:numel(elt(V1(end):end))-1
      single{1}(3+l) = elt(V1(end)-l);
    end
    
    if abs(V1(1)-V1(end)) > 1        % Verify if indices form an element
      single{2} = elt(V1(1):V1(end));
    end
  %end
end
