function PlotMesh(mesh, solution)
% AUTEUR : Diallo Amadou, 28/09/2020
%title('Approximate Solution');
max_n_vertices = max(cellfun(@length, mesh.cell_e));
padding_function = @(vertex_list) [vertex_list ...
			NaN(1,max_n_vertices-length(vertex_list))];
cell_e = cellfun(padding_function, mesh.cell_e, 'UniformOutput', false);
cell_e = vertcat(cell_e{:});
%data = [mesh.vertex, solution];
patch('Faces',cell_e,'Vertices',mesh.vertex,'FaceVertexCData', solution,'FaceColor','flat');
axis('square')
xlim([min(mesh.vertex(:,1)) - 0.1, max(mesh.vertex(:,1)) + 0.1])
ylim([min(mesh.vertex(:,2)) - 0.1, max(mesh.vertex(:,2)) + 0.1])
%zlim([min(solution) - 0.1, max(solution) + 0.1])
xlabel('x'); ylabel('y'); zlabel('u');
colormap('jet');
colorbar
end