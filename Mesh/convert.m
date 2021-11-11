
% convert to txt
%clear all
%close all
addpath('mesh_files/');

filename1 = "mesh1";
load(strcat(filename1, ".mat"));
filelocation = "../benchmark/";
fid = fopen(strcat(filelocation, filename, ".msh"), "w");

fprintf(fid, '$MeshFormat\n');
fprintf(fid, '2.2  0  8\n');
fprintf(fid, '$EndMeshFormat\n');
% vertices
fprintf(fid, '$Nodes\n');
Nvtx = size(vertices, 1);
fprintf(fid, '%d\n', Nvtx);
vertices = vertices';
for i = 1:Nvtx
    fprintf(fid, '%d %32.16f %32.16f\n', i, vertices(:,i));
end
fprintf(fid, '$EndNodes\n');
tasks = nbT*ones(length(elements),1);
lelem = fix(length(elements)/nbT);
for k = 1:nbT
    lelem_first = (k-1)*lelem+1;
    lelem_end = k*lelem;
    for i = lelem_first:lelem_end
        tasks(i) = k;
    end
end
% elements
fprintf(fid, '$Elements\n');
fprintf(fid, '%d\n', length(elements));
for i=1:length(elements)
    ne = length(elements{i});
    %fprintf(fid, ['%d %d %d 1 1 1 1 ', repmat('%d ', 1, ne), '\n'], i, ne, tasks(i), elements{i});
    fprintf(fid, ['%d %d  %d  ', repmat('%d ', 1, ne), '\n'], i, ne, tasks(i), elements{i});
end
fprintf(fid, '$EndElements\n');

% element to edges : e2edg
fprintf(fid, '$E2edg\n');
fprintf(fid, '%d\n', length(e2edg));
for i=1:length(e2edg)
    nedg = length(e2edg{i});
    %fprintf(fid, ['%d %d %d 1 1 1 1 ', repmat('%d ', 1, nedg), '\n'], i, nedg, tasks(i), e2edg{i});
    fprintf(fid, ['%d %d  %d  ', repmat('%d ', 1, nedg), '\n'], i, nedg, tasks(i), e2edg{i});
end
fprintf(fid, '$EndE2edg\n');

% neighbord of each element : e2neigh
fprintf(fid, '$E2neigh\n');
fprintf(fid, '%d\n', length(e2neigh));
for i=1:length(e2neigh)
    nneigh = length(e2neigh{i});
    %fprintf(fid, ['%d %d %d 1 1 1 1 ', repmat('%d ', 1, nneigh), '\n'], i, nneigh, tasks(i), e2neigh{i});
    fprintf(fid, ['%d %d  %d  ', repmat('%d ', 1, nneigh), '\n'], i, nneigh, tasks(i), e2neigh{i});
end
fprintf(fid, '$EndE2neigh\n');

% connectics (normal components) of each edges  : connect
fprintf(fid, '$Connect\n');
fprintf(fid, '%d\n', length(connect));
for i=1:length(connect)
    ncon = length(connect{i});
    %fprintf(fid, ['%d %d %d 1 1 1 1 ', repmat('%d ', 1, ncon), '\n'], i, ncon, tasks(i), connect{i});
    fprintf(fid, ['%d %d  %d  ', repmat('%d ', 1, ncon), '\n'], i, ncon, tasks(i), connect{i});
end
fprintf(fid, '$EndConnect\n');

% edges from the boundary : freeedg
fprintf(fid, '$Freeedg\n');
fprintf(fid, '%d\n', length(freeedg));
fprintf(fid, '%d\n', freeedg);
fprintf(fid, '$EndFreeedg\n');

%boundary
fprintf(fid, '$Boundary\n');
fprintf(fid, '%d\n', length(boundary));
fprintf(fid, '%d\n', boundary);
fprintf(fid, '$EndBoundary\n');

%Interior nodes
fprintf(fid, '$IntNodes\n');
fprintf(fid, '%d\n', length(freenod));
fprintf(fid, '%d\n', freenod);
fprintf(fid, '$EndIntNodes\n');

fclose(fid);