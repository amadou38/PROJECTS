function plot_solution(mesh, solution)
% AUTEUR : Diallo Amadou, 28/09/2020
title('Generated mesh');
max_n_vertices = max(cellfun(@length, mesh.elements));
padding_function = @(vertex_list) [vertex_list' ...
			NaN(1,max_n_vertices-length(vertex_list))];
elements = cellfun(padding_function, mesh.elements, 'UniformOutput', false);
elements = vertcat(elements{:});
%data = [mesh.vertices, solution];
%patch('Faces',elements,'Vertices',mesh.vertices,'FaceVertexCData', solution,'FaceColor','flat');
patch('Faces',elements,'Vertices',mesh.vertices,'FaceColor','none');
axis('square')
%grid on
title(['Generated mesh. Nelem: ', num2str(numel(mesh.elements)), ' elems & ', num2str(numel(mesh.vertices(:,1))), ' nodes']);
xlim([min(mesh.vertices(:,1)) - 0.1, max(mesh.vertices(:,1)) + 0.1])
ylim([min(mesh.vertices(:,2)) - 0.1, max(mesh.vertices(:,2)) + 0.1])
%zlim([min(solution) - 0.1, max(solution) + 0.1])
xlabel('x'); ylabel('y'); zlabel('u');
colormap('jet');
%colorbar
end