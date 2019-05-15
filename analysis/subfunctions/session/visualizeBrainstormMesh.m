function visualizeBrainstormMesh(anatDir, data, cmap, thresh, clims, meshType, ttl)



%% Check inputs

if ~exist('meshType','var') || isempty(meshType)
    bs_pial_low = load(fullfile(anatDir, 'tess_cortex_pial_low.mat'));
else
    if strcmp(meshType,'smooth')
        bs_pial_low = load(fullfile(anatDir, 'tess_cortex_pial_low_fig.mat'));
    elseif strcmp(meshType,'unsmooth')
        bs_pial_low = load(fullfile(anatDir, 'tess_cortex_pial_low.mat'));
    end
end


if ~exist('cmap','var') || isempty(cmap)
    cmap = hsv(256);
    
end


if ~exist('thresh','var') || isempty(thresh)
    thresh = 0;
end

if ~exist('clims','var') || isempty(clims)
    clims = [min(data(:)), max(data(:))];
end

if ~exist('ttl','var') || isempty(ttl)
    ttl = 'Brainstorm Mesh';
end

% Plot mean amplitude across epochs
figure; set(gcf, 'Color', 'w', 'Position', [163 483 891 554])

% Set up mesh
tH = trimesh(bs_pial_low.Faces,bs_pial_low.Vertices(:,1),bs_pial_low.Vertices(:,2),bs_pial_low.Vertices(:,3));
axis equal; hold on

% Preallocate space for colors to plot (nr vertices x 3 for RGB)
colors = NaN(size(bs_pial_low.Curvature,1),3);
sz     = length(cmap)-1;

% Get curvature
curv = bs_pial_low.Curvature; 

% Implement colors in curvature
colors(curv<=0,:) = .25;
colors(curv>0,:) = .75;

% Get index for data above the requested thresh (default = 0) and select
% those data
ii = find(data>thresh);
Z = data(ii);

% Convert to 1-266
Z_ind = round(sz.*((Z-min(Z)) ./ (max(Z)-min(Z))))+1;

% overlay in colors variable
colors(ii,:) = cmap(Z_ind,:);

% set source pediction as colors
set(tH, 'LineStyle', 'none', 'FaceColor', 'interp', 'FaceVertexCData', colors);
% 
colormap(cmap); colorbar; set(gca, 'CLim',clims);

pos = [-.1 0 .1];
light('Position',pos,'Style','local')
% lighting gouraud
material shiny; %dull
title(sprintf('%s', ttl)); 

drawnow;


return