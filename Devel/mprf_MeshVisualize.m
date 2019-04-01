% (1) Build a freesurfer compatible mesh on mrMesh
% (2) Load the saved prf parameter values to that

% building mrMesh

FS_surface = '/mnt/storage_2/MEG/Retinotopy/Data/Freesurfer_directory/wlsubj004/surf';
fs_prf_data = '/mnt/storage_2/MEG/Retinotopy/Quality_check/wlsubj004/prf_data/surface/freesurfer';

surfaces_to_load = {'lh.pial','rh.pial'};

n=1;
cur_surf = fullfile(FS_surface,surfaces_to_load{n});

tmp = strsplit(surfaces_to_load{n},'.');
cur_hs = tmp{1};
fprintf('Exporting parameters for %s hemisphere:\n',cur_hs);

% Load mesh using fs_meshFromSurface, this creates a mrVista compatible
% mesh. Using 'my own' function that skips the smoothing:
mrv_msh = mprf_fs_meshFromSurface(cur_surf);
fnum = numel(mrv_msh.triangles);

% compute mapping using mrmMapVerticesToGray (mrmMapGrayToVertices):
cur_v2gmap = mrmMapVerticesToGray(mrv_msh.vertices, viewGet(hvol,'nodes'),...
    mmPerVox);

% Add the vertexGrayMap field to mesh properties
mrv_msh = meshSet(mrv_msh,'vertexGrayMap',cur_v2gmap);


%% Making mesh from freesurfer surface ELine's code
mshFS = fs_meshFromSurface(cur_surf);
mshFS = meshSmooth(mshFS);
mshFS = meshColor(mshFS);

data_type = 'Averages';
tmp_sm_par = curvWangAtlas;

hvol = initHiddenGray;
hvol = viewSet(hvol,'curdt',data_type);
hvol = viewSet(hvol,'map',{tmp_sm_par'});


mmPerVox = viewGet(hvol,'mmpervox');
% compute mapping using mrmMapVerticesToGray (mrmMapGrayToVertices):
cur_v2gmap = mrmMapVerticesToGray(mshFS.vertices, viewGet(hvol,'nodes'),...
    mmPerVox);

% Add the vertexGrayMap field to mesh properties
mshFS = meshSet(mshFS,'vertexGrayMap',cur_v2gmap);
% Add mesh to the volume view:
hvol = viewSet(hvol,'add and select mesh',mshFS);

% meshVisualize(mshFS);

% Update mesh
%? hvol = meshColorOverlay(hvol);

subjectSessionDir = '/mnt/storage_2/MEG/Retinotopy/Subject_sessions';
subject = 'wlsubj004';
hemi = 'lh';
wangAtlasPth = fullfile(subjectSessionDir, subject, 'rois', 'surface', 'freesurfer', sprintf('%s.V1', hemi));
[curvWangAtlas, fnum] = read_curv(wangAtlasPth);
 
figure(1)    
% Faces (also called triangles) are defined by 3 points, each of
% which is an index into the x, y, z vertices
faces = mshFS.triangles' + 1; % we need to 1-index rather than 0-index for Matlab
    
% The vertices are the locations in mm spacing
x     = mshFS.vertices(1,:)';
y     = mshFS.vertices(2,:)';
z     = mshFS.vertices(3,:)';
    
% The colormap will, by default, paint sulci dark and gyri light
c     = mshFS.colors(1,:)';
    
% Render the triangle mesh
tH = trimesh(faces, x,y,z);
 
% Make it look nice
set(tH, 'LineStyle', 'none', 'FaceColor', 'interp', 'FaceVertexCData',c)
axis equal off; colormap gray; set(gca, 'CLim', [0 255])
    
% Lighting to make it look glossy
light('Position',100*[0 1 1],'Style','local')
% lighting gouraud
 
% Which mesh are we plotting?
title(sprintf('Freesurfer mesh: %s', hemi))
 
% Rotate it
set(gca, 'View', [-16.7000  -90.0000]);


% find the vertices that belong to the ROI 
% Determine the faces associated with it






hold on;



% Render the triangle mesh
tH_cm = trimesh(faces, x(roiIdx),y(roiIdx),z(roiIdx));

% The colormap will, by default, paint sulci dark and gyri light
cm     = mshFS.colors(1,:)';
nanIdx = isnan(curvWangAtlas);
roiIdx = find(curvWangAtlas(~nanIdx));
cm(roiIdx) = curvWangAtlas(~nanIdx);
set(tH_cm, 'LineStyle', 'none', 'FaceColor', 'interp', 'FaceVertexCData',cm)
axis equal off; colormap default; set(gca, 'CLim', [0 255])









%%

good_mapping = cur_v2gmap > 0;




% Build mesh
hvol = meshBuild(hvol,'left');

% 
MSH = meshVisualize(viewGet(hvol,'Mesh')); hvol = viewSet(hvol, 'Mesh', MSH); clear MSH;

% Smooth the mesh
hvol = viewSet(hvol, 'Mesh', meshSmooth( viewGet(hvol, 'Mesh'), 1));

% Load Roi local. for shared change the second 1 to 0
hvol = loadROI(hvol, 'dialog', [],[],1,1); hvol = selectCurROISlice(hvol); hvol = refreshScreen(hvol,0);

% Update mesh
hvol = meshColorOverlay(hvol);



%%
data_type = 'Averages';
hvol = initHiddenGray;
hvol = viewSet(hvol,'curdt',data_type);

% Build mesh
hvol = meshBuild(hvol,'left');

% Load the beta map which was projected on to freesurfer space
tmp = read_curv('/mnt/storage_2/MEG/Retinotopy/Subject_sessions/wlsubj004/prf_data/surface/freesurfer/lh.beta');
idx_nan = isnan(tmp);
tmp(idx_nan) = 0;

map = tmp;

hvol = viewSet(hvol,'map',{tmp'});

hvol = meshColorOverlay(hvol);


