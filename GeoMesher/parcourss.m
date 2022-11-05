%function pelt = parcourss(elt, S)
clear; close all;
addpath('Connectique/');
addpath('mesh_files/');
num_ind = 6;   % Corresponding face indices in the 'ind' cell array. NB: good from 1 to 13 
N = 1;        % N^2 corresponds to the total of principal elements on the mesh       
nbT = 1;
L = [1 1];     % DOM = [0 L(1)]x[0 L(2)]
ind = indice2();    % function containing the default indices of principal elements
N = N*(max(ind{num_ind}(1,:))-1); % size of the box that contain elements 

h1 = (L(1)-0)/N;
h2 = (L(2)-0)/N;
[X,Y] = meshgrid(0:h1:L(1),0:h2:L(2));
vertices = [X(:), Y(:)];
[n,m] = size(X);
A = reshape(1:n*m, n, m);
nx = n; ny = m;
ind = ind{num_ind}
ni = max(ind(1,:)); mj = max(ind(2,:)); 
ni = ni - 1; mj = mj - 1;
n = fix((n-1)/ni); 
m = fix((m-1)/mj);
% First principal element
x = 1:max(ind(1,:)); y = 1:max(ind(2,:));
S = A(x,y);
elt = zeros(1,numel(ind(1,:)));
for j = 1:numel(ind(1,:))
    elt(j) = S(ind(1,j), ind(2,j));
end

  %for i = 1:numel(elt)
    
    s = cell(4,1); sing = s; Vi = s; Ii = s; Belt = cell(1,1);
    s{1} = S(:,1); s{2} = S(end,:); s{3} = S(:,end); s{4} = S(1,:);
    count = 0;
    for i = 1:4
        sing{i}(1) = s{i}(1);
        [Ii{i},Vi{i}] = find(ismember(elt,s{i}));
        k = 2;
        if elt(Vi{i}(1)) ~= s{i}(1)
            k = 1;
        end
        if numel(Ii{i}) > 1
            sing{i}(2:numel(Vi{i}(k:end))+1) = elt(Vi{i}(k:end));
        end
        if elt(Vi{i}(end)) ~= s{i}(end)
            sing{i}(end+1) = s{i}(end);
        end
        if (i >= 3)
            sing{i} = sing{i}(end:-1:1);
        end
        sing{i} = unique(sing{i},'stable');
        Belt{1}(1) = sing{i}(1); p = 2;
        if sing{i}(1) == elt(Vi{i}(1))
            p = 1;
        end
        for id = Vi{i}(1):Vi{i}(end)
            Belt{1}(p) = elt(id);
            p = p + 1;
        end
        if sing{i}(end) ~= elt(Vi{i}(end))
            Belt{1}(p) = sing{i}(end);
        end
        Belt{1} = unique(Belt{1},'stable'); 
        Belt{1}
        count = count + abs(numel(Belt{1}) - numel(sing{i}));    
    end
    
    
    
    
    single{1}(1) = S(1,1);
    single = cell(1,1);  
    [I,V1] = find(ismember(elt,S(:,1)));
    if numel(I) > 1
      single{1}(2:numel(V1(2:end))+1) = elt(V1(2:end));
    end
    sing{:}
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
  
elt

%   single1 = left(elt, S);
%   single2 = top(elt, S);
%   single3 = right(elt, S);
%   single4 = bottom(elt, S);
%   %pelt = cell(1, numel(single1)+numel(single2)+numel(single3)+numel(single4));
%   pelt = {single1{1}, single2{1}, single3{1}, single4{1}, single1{2:end}, single2{2:end}, single3{2:end}, single4{2:end}};
% %end
