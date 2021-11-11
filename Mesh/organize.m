function [elts,vertices,boundary] = organize(elts, vertices)

nod = unique(cell2mat(elts));
vect = 1:numel(nod);
elts1 = elts;
for k = 1:length(nod)
    for i = 1:numel(elts)
      for j = 1:numel(elts{i})
        if elts{i}(j) == nod(k)
          elts1{i}(j) = vect(k);
        end
      end
    end
end

elts = elts1;
vertices = vertices(ismember(1:size(vertices,1),nod),:);

boundary = unique([find(ismember(vertices(:,1),min(vertices(:,1)))); find(ismember(vertices(:,2),max(vertices(:,2)))); find(ismember(vertices(:,1),max(vertices(:,1)))); find(ismember(vertices(:,2),min(vertices(:,2))))]);


%connect = alg6(e2edg,nedg,nelt);

end