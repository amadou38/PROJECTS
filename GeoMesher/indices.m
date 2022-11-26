function ind = indices(nb_pelt)
n = 1;
indices{n} = [1 2 1; 2 2 1]; n = n+1;
indices{n} = [1 2 2; 2 2 1]; n = n+1;

indices{n} = [1 1 2 1 1 2 3 4 5 5 4 5 5 4 3 2; 1 2 3 4 5 5 4 5 5 4 3 2 1 1 2 1]; n = n+1;
indices{n} = [1 2 1 3 5 4 5 3; 1 3 5 4 5 3 1 2]; n = n+1;
indices{n} = [1 2 1 2 3 4 5 4 5 3; 1 3 4 4 5 4 4 3 1 2]; n = n+1; % stars
indices{n} = [1 1 2 3 3 2; 1 2 3 3 2 1];                          % hexagon

indices{n} = [1 2 2 1; 2 2 1 1]; n = n+1;
indices{n} = [1 2 3 2; 2 3 2 1]; n = n+1;
indices{n} = [1 2 2 1 2 3 4 5 6 5 5 6 5 4 3 2; 2 3 4 5 6 5 5 6 5 4 3 2 1 2 2 1]; n = n+1;
indices{n} = [1 3 5 4 2; 3 5 3 1 1]; n = n+1;                    % pentagon
indices{n} = [1 2 3 4 5 4 2; 2 4 5 4 2 1 1]; n = n+1;              % heptagon
indices{n} = [1 1 2 3 4 4 3 2; 2 3 4 4 3 2 1 1]; n = n+1;        % octogon
indices{n} = [1 3 1 2 4 5 5 4 2; 2 3 4 5 5 4 2 1 1]; n = n+1;
indices{n} = [1 1 3 5 5 3; 2 4 5 4 2 1]; n = n+1;
indices{n} = [1 3 4 5 7 5 4 3; 4 5 7 5 4 3 1 3]; n = n+1;
indices{n} = [1 2 3 2; 2 3 2 1]; n = n+1;
indices{n} = [1 3 3 5 3; 3 5 3 3 1]; n = n+1;
indices{n} = [1 1 2 3 4 2 2 4 3 2; 2 3 4 4 3 3 2 2 1 1]; n = n+1;
indices{n} = [1 2 4 5 4 2; 3 5 5 3 1 1]; n = n+1;
indices{n} = [1 2 1 2 3 4 5 4 5 4 3 2; 2 3 4 5 4 5 4 3 2 1 2 1]; n = n+1; % X-form for test
indices{n} = [1 2 3 2 1 2 3 4 5 4 4 5 4 3 2; 2 2 3 4 4 5 4 5 4 4 2 2 1 2 1]; n = n+1; % X-form for test2


ind = indices{nb_pelt};
ind = ind(end:-1:1,:);
end
