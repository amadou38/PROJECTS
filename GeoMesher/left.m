function single = left(elt,S)
  
  %for i = 1:numel(elt)
    single = cell(1,1);
    if ~ismember(S(1,1),elt)
      single{1}(1) = S(1,1);
    end
    [I,V1] = find(ismember(elt,S(:,1)));
    if I ~= 0
      single{1}(2) = elt(V1(1));
    end
    [I,V2] = find(ismember(elt,S(1,:)));
    for l = 1:numel(elt(V2(end):end))-1
      single{1}(2+l) = elt(V2(1)+l);
    end
    if I ~= 0
      single{1}(end+1) = elt(V2(end));
    end
    
    if abs(V1(1)-V1(end)) > 1        % Verify if indices form an element
      single{2} = elt(V1(1):V1(end));
    end
  %end
%   if rep == 1
%     this = single{1};
%     single{1}(end) = this(1);
%     single{1}(1:end-1) = this(2:end);
%   end  
end
