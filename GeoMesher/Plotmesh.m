function Plotmesh(mesh,U)

vertices = mesh.vertices;
elements = zeros(length(U),size(mesh.elements{1},2));
for i = 1:length(U)
    elements(i,:) = mesh.elements{i};
end

patch('Faces',elements,'Vertices',vertices,'FaceVertexCData', U,'FaceColor','flat');
axis('square')
xlabel('x'); ylabel('y'); zlabel('u');
colormap('jet');

end