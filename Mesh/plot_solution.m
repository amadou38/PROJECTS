function plot_solution(mesh,tit, solution)
% AUTEUR : Diallo Amadou, 28/09/2020
title(tit);
max_n_vertices = max(cellfun(@length, mesh.elements));
padding_function = @(vertex_list) [vertex_list' ...
			NaN(1,max_n_vertices-length(vertex_list))];
elements = cellfun(padding_function, mesh.elements, 'UniformOutput', false);
elements = vertcat(elements{:});
data = [mesh.vertices, solution];
for i = 1:length(solution)
    if data(i,3) == 3
        solution(i) = 0;
        data(i,3) = 'w';
    end
end
patch('Faces', elements,...
	'Vertices', data,...
	'FaceColor', 'interp',... 
	'CData', solution);
patch('Faces',elements,'Vertices',mesh.vertices,'FaceColor','none');
axis('square')
xlim([min(mesh.vertices(:,1)) - 0.1, max(mesh.vertices(:,1)) + 0.1])
ylim([min(mesh.vertices(:,2)) - 0.1, max(mesh.vertices(:,2)) + 0.1])
zlim([min(solution) - 0.1, max(solution) + 0.1])
xlabel('x'); ylabel('y'); zlabel('u');
colormap('jet');
colorbar
end