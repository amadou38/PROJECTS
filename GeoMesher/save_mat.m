function save_mat(name, elements, nod2e, vertices, e2edg, e2neigh, connect, freeedg, boundary,freenod)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

filename = ['mesh_files/', name, '.mat'];
save(filename, 'elements', 'nod2e', 'vertices', 'e2edg', 'e2neigh', 'connect', 'freeedg', 'boundary','freenod');


end

