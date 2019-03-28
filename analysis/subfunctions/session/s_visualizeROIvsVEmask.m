% s_visualizeROIvsVEmask

subjectSessionDir = '/Volumes/server/Projects/MEG/Retinotopy/Subject_sessions';
bsDir = '/Volumes/server/Projects/MEG/brainstorm_db';
subject = 'wlsubj058';
anatDir = fullfile(bsDir, 'MEG_Retinopy', 'anat', subject);

cd(fullfile(subjectSessionDir, subject, 'mask'))


tmp = load('roimask.mat');
m.roi = double(tmp.cur_in);
tmp = load('vemask_withroimask.mat');
m.both = double(tmp.cur_in);
tmp = load('vemask_noroimask.mat');
m.ve = double(tmp.cur_in);

clear tmp;



visualizeBrainstormMesh(anatDir, m.both, [],[],[], 'ROI+VE mask');
print(gcf,'-dpng', [] ,fullfile(subjectSessionDir, subject, 'mask', 'bs_mesh_roive'));
savefig(gcf,fullfile(subjectSessionDir, subject, 'mask', 'bs_mesh_roive.fig'));

visualizeBrainstormMesh(anatDir, m.ve, [],[],[], 'VE mask');
print(gcf,'-dpng', [] ,fullfile(subjectSessionDir, subject, 'mask', 'bs_mesh_ve'));
savefig(gcf,fullfile(subjectSessionDir, subject, 'mask', 'bs_mesh_ve.fig'));

visualizeBrainstormMesh(anatDir, m.roi, [],[],[], 'ROI mask');
print(gcf,'-dpng', [] ,fullfile(subjectSessionDir, subject, 'mask', 'bs_mesh_roi'));
savefig(gcf,fullfile(subjectSessionDir, subject, 'mask', 'bs_mesh_roi.fig'));


fprintf('Number of total vertices ROI+VE mask: %d\n', sum(m.both));
fprintf('Number of total vertices VE mask: %d\n', sum(m.ve));
fprintf('Number of total vertices ROI mask: %d\n', sum(m.roi));
fprintf('Number of vertices overlapping ROI+VE vs VE only mask: %d\n', length(intersect(find(m.both),find(m.ve))));
fprintf('Number of vertices overlapping ROI+VE vs ROI only mask: %d\n', length(intersect(find(m.both),find(m.roi))));
fprintf('Number of vertices overlapping ROI only vs VE only mask: %d\n', length(intersect(find(m.roi),find(m.ve))));


%% Load Freesurfer meshes

% see t_initAnatomyFromFreesurfer, meshImportFreesurferSurfaces
hemi = 'lh';
fsHemiPth    = fullfile('/Volumes/server/Freesurfer_subjects/', subject, 'surf', sprintf('%s.inflated',hemi));

% Load mesh from FreeSurfer
mshFS = fs_meshFromSurface(fsHemiPth);


% Load data from Wang Atlas
wangAtlasPth = fullfile(subjectSessionDir, subject, 'rois', 'surface', 'freesurfer', sprintf('%s.all_rois', hemi));
[curvWangAtlas, fnum] = read_curv(wangAtlasPth);

% Visualize with vistasoft
% meshVisualize(mshFS)



% Visualize with mMatlab
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

nanIdx = isnan(curvWangAtlas);
roiIdx = find(curvWangAtlas(~nanIdx));
c(roiIdx) = 255;
    
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
