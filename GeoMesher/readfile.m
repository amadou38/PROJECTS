
% convert to txt
%clear all
%close all
addpath('mesh_files/');

load(strcat(filename, ".msh"));
% filelocation = "../benchmark/";
fid0 = fopen(strcat(filename, ".msh"), "r");
fid1 = fopen(strcat(filename, ".msh"), "r");

tasks = nbT*ones(length(elements),1);
lelem = fix(length(elements)/nbT);
for k = 1:nbT
    lelem_first = (k-1)*lelem+1;
    lelem_end = k*lelem;
    for i = lelem_first:lelem_end
        tasks(i) = k;
    end
end
Tasks = zeros(length(vertices(:,1)),1);
for i = 1:length(elements)
    for j = 1:length(elements{i})
        Tasks(elements{i}(j)) = tasks(i);
    end
end
formatSpec = '%f';
u0 = fscanf(fid0, formatSpec);
u1 = fscanf(fid1, formatSpec);
U = 3*ones(length(vertices(:,1)),1);
k = 1; l = 1;
for i = 1:length(U)
    if Tasks(i) == 1 && k <= length(u0)
        U(i) = u0(k);
        k = k + 1;
    elseif Tasks(i) == 2 && l <= length(u1)
        U(i) = u1(l);
        l = l + 1;
    end
end

fclose(fid0);
fclose(fid1);