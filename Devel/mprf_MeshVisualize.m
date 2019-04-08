%% mprf_MeshVisualize.m
%
% Example script to visualize rois and saved prf parameters on subject
% mesh, either with mrVista, on the Brainstorm mesh, or on the FS mesh
%
% (1) Visualize Wang atlas on FreeSurfer surface with mrVista
% (2) Visualize Wang atlas on BrainStorm surface
% (3) Visualize exported prf data on FreeSurfer surface without mrVista
% (4) Visualize exported prf data on Brainstorm surface


%% (0) General parameters
subject = 'wlsubj004';

fsDir = sprintf('/Volumes/server/Freesurfer_subjects/%s/', subject);
sessionDir = sprintf('/Volumes/server/Projects/MEG/Retinotopy/Subject_sessions/%s/', subject);
vistaSessionDir = sprintf('/Volumes/server/Projects/MEG/Retinotopy/Data/fMRI/%s/vistaSession/', subject);
brainstormAnatDir = sprintf('/Volumes/server/Projects/MEG/brainstorm_db/MEG_Retinotopy/anat/%s', subject);


roiDataDir = fullfile(sessionDir, 'rois', 'surface');
prfDataDir = fullfile(sessionDir, 'prf_data', 'surface');

surfaces_to_load = 'pial';
hemi = {'lh','rh'};

%% (1) Visualize Wang atlas on FS surface in mrVista

fprintf('Loading Wang atlas ROIs on mrVista hemispheres\n');

cd(vistaSessionDir);

% Load meshes from Vista
vw = mrVista('3');
mesh1 = fullfile('3DAnatomy', 'Left', '3DMeshes', 'Left_inflated.mat');
mesh2 = fullfile('3DAnatomy', 'Right', '3DMeshes', 'Right_inflated.mat');

[vw, OK] = meshLoad(vw, mesh1, 1); if ~OK, error('Mesh server failure'); end
[vw, OK] = meshLoad(vw, mesh2, 1); if ~OK, error('Mesh server failure'); end

% Define atlas path
wangAtlasPthNifti = fullfile(fsDir, 'mri', sprintf('native.wang2015_atlas.nii.gz'));
    
% Convert mgz to nifti   
vw = nifti2ROI(vw, wangAtlasPthNifti);

% Let's look at the ROIs on meshes
vw = roiSetVertIndsAllMeshes(vw); 

nROIs = length(viewGet(vw, 'ROIs'));
colors = hsv(nROIs);
for ii = 1:nROIs
   vw = viewSet(vw, 'ROI color', colors(ii,:), ii); 
end
vw = viewSet(vw, 'roi draw method', 'boxes');

vw = meshUpdateAll(vw); 


%% (2) Visualize Wang atlas on BS surface

fprintf('Loading Wang atlas ROIs on Brainstorm mesh:\n');

wangAtlasPthBS = fullfile(roiDataDir, 'brainstorm', 'pial.all_rois');
[curvWangAtlas, fnum] = read_curv(wangAtlasPthBS);

ttl = 'Brainstorm Wang atlas ROIs';
visualizeBrainstormMesh(brainstormAnatDir, curvWangAtlas, [], [], [], ttl)
colorbar;

%% (3) Visualize prf data on FS surface

prfParam = 'beta';

for h = 1:length(hemi)
    fprintf('Loading prf data: %s on FreeSurfer mesh: %s\n', hemi{h}, prfParam);

    
    fsHemiPth    = fullfile(fsDir, 'surf', sprintf('%s.inflated',hemi{h}));
    
    % Load mesh from FreeSurfer
    mshFS = fs_meshFromSurface(fsHemiPth);

    % Load data from Wang Atlas
    prfDataPth = fullfile(prfDataDir, 'freesurfer', sprintf('%s.%s', hemi{h}, prfParam));

    [curvPrfData, fnum] = read_curv(prfDataPth);

    % Visualize with Matlab
    figure;    
    
    % Faces (also called triangles) are defined by 3 points, each of
    % which is an index into the x, y, z vertices
    faces = mshFS.triangles' + 1; % we need to 1-index rather than 0-index for Matlab

    % The vertices are the locations in mm spacing
    x     = mshFS.vertices(1,:)';
    y     = mshFS.vertices(2,:)';
    z     = mshFS.vertices(3,:)';

    % The colormap will, by default, paint sulci dark and gyri light
    cmap = [gray(128); jet(128)];
    curvature = mshFS.colors(1,:)';
    colors = curvature;
    colors(mshFS.curvature<0) = -1.5;
    colors(mshFS.curvature>=0) = -.5;
    
    nanIdx = isnan(curvPrfData);

    dataIdx = find(curvPrfData(~nanIdx));

    colors(dataIdx) = curvPrfData(dataIdx);
%     colors(nanIdx) = -.5;
    % Render the triangle mesh
    tH = trimesh(faces, x,y,z);
    
    % Make it look nice
    set(tH, 'LineStyle', 'none', 'FaceColor', 'interp', 'FaceVertexCData',colors)
    axis equal off; colormap(cmap); set(gca, 'CLim', [-2, 2])
    colorbar;
    
    % Lighting to make it look glossy
    light('Position',100*[0 1 1],'Style','local')

%     lighting gouraud

    % Which mesh are we plotting?
    title(sprintf('Freesurfer mesh: %s', hemi{h}))

    % Rotate it
    set(gca, 'View', [-16.7000  -90.0000]);

end

%% (4) Visualize prf data on BS surface

prfParam = 'beta';

fprintf('Loading prf data: %s on Brainstorm mesh: %s\n', hemi{h}, prfParam);


prfDataPthBS = fullfile(prfDataDir, 'brainstorm', sprintf('pial.%s', prfParam));
[curvPrfBS, fnum] = read_curv(prfDataPthBS);

nanIdx  = isnan(curvPrfBS);
% dataIdx = find(curvPrfBS);

curvPrfBS(nanIdx) = -0.5;

colors = curvPrfBS;
clims = [-1 1];

thresh = 0;

ttl = sprintf('Brainstorm prf data: %s', prfParam);
visualizeBrainstormMesh(brainstormAnatDir, colors, thresh, clims, [], ttl)

colorbar;


