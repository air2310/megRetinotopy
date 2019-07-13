function bs_msh = mprfMeshFromBrainstorm(surface_file)

bs_surf = load(surface_file);

vertices = bs_surf.Vertices' .* 1000;
faces = bs_surf.Faces(:,[2 1 3])';

bs_msh = meshCreate;

bs_msh.initVertices = vertices;
bs_msh.triangles = faces-1;   % vista triangles are zero indexed
bs_msh.vertices = bs_msh.initVertices;
bs_msh.colors = ones(4, size(bs_msh.vertices,2))*128;

% Construct surface normals
% TR = triangulation(bs_msh.triangles'+1, bs_msh.vertices');
% VN = vertexNormal(TR);
bs_msh.normals = bs_surf.VertNormals';


bs_msh.smooth_iterations = 15;
bs_msh.smooth_relaxation = 1;
bs_msh = meshSmooth(bs_msh);
bs_msh = meshColor(bs_msh);

end








