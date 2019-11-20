%% Plot dipoles on mesh

% Define subject
subjID   = 'wlsubj111'; 

% Load paths with data files for this subject
dirPth = loadPaths(subjID);

% Go back to root
cd(mprf_rootPath)

% Set options
opt = getOpts('saveFig',1,'verbose',1); % see getOpts function for more options


% Brainstorm pial mesh
bs_pialSurface = load(fullfile(dirPth.bs.anatPth, 'tess_cortex_pial_low.mat'));

% Freesurfer pial mesh
fs_pialSurface = load(fullfile(dirPth.bs.anatPth, 'tess_cortex_pial_high.mat'));

% colors = abs(template.V1StimEccen);
    
% Define colorbar colors
cmap = [gray(128)];

%% Plot Brainstorm mesh
figure; set(gcf, 'Color', 'w', 'Position', [163 483 891 554])

% Set up Brainstorm mesh
tH = trisurf(bs_pialSurface.Faces,bs_pialSurface.Vertices(:,1),bs_pialSurface.Vertices(:,2),bs_pialSurface.Vertices(:,3));
axis equal; hold on

quiver3(bs_pialSurface.Vertices(:,1),bs_pialSurface.Vertices(:,2),bs_pialSurface.Vertices(:,3), ...
bs_pialSurface.VertNormals(:,1),bs_pialSurface.VertNormals(:,2),bs_pialSurface.VertNormals(:,3))

% set source pediction as colors
set(tH, 'LineStyle', 'none', 'FaceColor', 'interp');%, 'FaceVertexCData',double(colors));

colormap(cmap); colorbar; set(gca, 'CLim', [-2 2]);
title(sprintf('Brainstorm pial mesh dipoles - %s', subjID));
pos = [-.1 0 .1];
light('Position',pos,'Style','local')
lighting gouraud
material shiny; %dull

%% Freesurfer mesh

% Plot mean amplitude across epochs
figure; set(gcf, 'Color', 'w', 'Position', [163 483 891 554])

% Set up Brainstorm mesh
tH = trisurf(fs_pialSurface.Faces,fs_pialSurface.Vertices(:,1),fs_pialSurface.Vertices(:,2),fs_pialSurface.Vertices(:,3));
axis equal; hold on

quiver3(fs_pialSurface.Vertices(:,1),fs_pialSurface.Vertices(:,2),fs_pialSurface.Vertices(:,3), ...
fs_pialSurface.VertNormals(:,1),fs_pialSurface.VertNormals(:,2),fs_pialSurface.VertNormals(:,3))

% set source pediction as colors
set(tH, 'LineStyle', 'none', 'FaceColor', 'interp');%, 'FaceVertexCData',double(colors));

colormap(cmap); colorbar; set(gca, 'CLim', [-2 2]);
title(sprintf('Freesurfer pial mesh dipoles - %s', subjID));
pos = [-.1 0 .1];
light('Position',pos,'Style','local')
lighting gouraud
material shiny; %dull