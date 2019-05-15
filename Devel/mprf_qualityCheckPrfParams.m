%% mprf_qualityCheckPrfParams.m
% Quality checks for multiple steps in the analysis pipeline:
%
% 1) Projection of pRF parameters and ROIs from each voxel (mrVista space)
%    to each vertex (freesurfer space)
% 2) ROI projection / Predicted responses from each vertex (stimulus * pRF models)
% 3) Projection of parameters / ROIs onto brainstorm space
% 4) Synthetic dataset to see if pipeline works (Measured signal =
%    Predicted signal + noise ?)
% 5) Plot all 19 reference phases for a sensor with high VE
% 6) VE maps and predicted responses for every size and position range

% By Eline Kupers, NYU 2019
%% 0. Define params and paths

subject       = 'wlsubj030';
projectFolder = '/Volumes/server/Projects/MEG/Retinotopy/';
saveFigFolder = fullfile(projectFolder, 'Quality_check', sprintf('%s',subject));
serverFolder  = fullfile(projectFolder, 'Data');
vistaDataFolder = fullfile(serverFolder, 'fMRI', subject, 'vistaSession');
vistaDataType   = 'Averages';
vistaRetModelFolder = fullfile(vistaDataFolder,'Gray', vistaDataType);
vistaClassFile  = fullfile(vistaDataFolder, '3DAnatomy', 't1_class.nii.gz');

d = dir(fullfile(vistaRetModelFolder, 'rm_*fFit.mat'));
vistaRetModel = d.name;

fsDir = fullfile('/Volumes/server/Freesurfer_subjects/',subject);
prfDataDir = fullfile('/Volumes/server/Projects/MEG/Retinotopy/Subject_sessions', subject, 'prf_data');
roiDataDir = fullfile('/Volumes/server/Projects/MEG/Retinotopy/Subject_sessions', subject, 'rois');

brainstormAnatDir = fullfile('/Volumes/server/Projects/MEG/brainstorm_db/MEG_Retinotopy/anat/',subject);

%% 1. Initial pRF parameters in mrVista space

% % Go to folder
% cd(vistaDataFolder);
% 
% % Load session and check anatomy path
% load('mrSESSION.mat')
% if strcmp(vANATOMYPATH, '3DAnatomy/vAnatomy.dat')
%     vANATOMYPATH = fullfile(vistaDataFolder, '3DAnatomy', 't1.nii.gz');
%     mrSESSION.vANATOMYPATH = vANATOMYPATH;
% end
% 
% % Load meshes from Vista
% vw = mrVista('3');
% mesh1 = fullfile('3DAnatomy', 'Left', '3DMeshes', 'Left_inflated.mat');
% mesh2 = fullfile('3DAnatomy', 'Right', '3DMeshes', 'Right_inflated.mat');
% 
% [vw, OK] = meshLoad(vw, mesh1, 1); if ~OK, error('Mesh server failure'); end
% [vw, OK] = meshLoad(vw, mesh2, 1); if ~OK, error('Mesh server failure'); end
% 
% % Set the volume view to the current data type and add the RM model
% vw = viewSet(vw,'curdt',vistaDataType);
% vw = rmSelect(vw,1,vistaRetModel);
% 
% % % Mask to exclude unreliable voxels (i.e. VE == 0) from the smoothing
% % % below. Otherwise, the pRF parameters will be averaged with a lot of
% % % zeros:
% % sm_mask = rmGet(hvol.rm.retinotopyModels{1},'varexplained') > 0;
% % % We need these parameters from the pRF model
% % params = {'sigma','x','y','varexplained','beta'};
% % % We need the mrVista segmentation to check if the selection of pRF
% % % parameters is correct, i.e. all selected pRF parameters must fall in the
% % % gray matter.
% % cls = niftiRead(vistaClassFile);
% 
% vw = meshUpdateAll(vw);

%% 2. Check original versus smoothed parameters (x,y,sigma, beta) on FS surface

prfParams = {'eccentricity', 'polar_angle', 'beta', 'eccentricity_smoothed','polar_angle_smoothed','recomp_beta', 'varexplained'};

hemi = {'lh', 'rh'};
% clim = [0 1]; % eccen, polar angle and beta color maps


%% Plot pRF params on FS surface

% load var expl separate to mask data
        
% Load data from surface folder
[varExpl.lh, ~] = read_curv(fullfile(prfDataDir, 'surface', 'freesurfer', 'lh.varexplained'));
[varExpl.rh, ~] = read_curv(fullfile(prfDataDir, 'surface', 'freesurfer', 'rh.varexplained'));


for p = 1:length(prfParams)
    
    for h = 1:length(hemi)
        fprintf('Loading prf data: %s on FreeSurfer mesh: %s\n', hemi{h}, prfParams{p});
        
        fsHemiPth    = fullfile(fsDir, 'surf', sprintf('%s.inflated',hemi{h}));
        
        % Load mesh from FreeSurfer
        mshFS = fs_meshFromSurface(fsHemiPth);
        
        % Load data from surface folder
        prfDataPth = fullfile(prfDataDir, 'surface', 'freesurfer', sprintf('%s.%s', hemi{h}, prfParams{p}));
        
        [curvPrfData, fnum] = read_curv(prfDataPth);
        nanIdx = isnan(curvPrfData);   
        figure; histogram(curvPrfData(~nanIdx)); ylabel('Frequency'); xlabel('vertex value');
        set(gca, 'TickDir', 'out', 'FontSize', 14)
        title(sprintf('Freesurfer mesh: %s %s', hemi{h}, prfParams{p}))
         if (strcmp(prfParams{p}, 'beta') || strcmp(prfParams{p} ,'recomp_beta'))
             mn = mean(curvPrfData(~nanIdx));
             sd3 = (3*std(curvPrfData(~nanIdx)));
             if isnan(sd3); sd3=0; end
            set(gca, 'XLim', [0, mn+sd3]);
        end
        
        print(gcf, '-dpng', fullfile(saveFigFolder,sprintf('%s_FSmesh_%s_%s_histogram', subject, hemi{h}, prfParams{p})))
        
        % Faces (also called triangles) are defined by 3 points, each of
        % which is an index into the x, y, z vertices
        faces = mshFS.triangles' + 1; % we need to 1-index rather than 0-index for Matlab
        
        % The vertices are the locations in mm spacing
        x     = mshFS.vertices(1,:)';
        y     = mshFS.vertices(2,:)';
        z     = mshFS.vertices(3,:)';
        
        % The colormap will, by default, paint sulci dark and gyri light
        cmap = hsv(266);
        
        % Get curvature
        curv = mshFS.curvature; 

        % Preallocate space for colors to plot (nr vertices x 3 for RGB)
        colors = NaN(size(curv,2),3);
        sz     = length(cmap)-1;
        clim   = [0 max(curvPrfData(:))];
        
        % Implement colors in curvature
        colors(curv<=0,:) = .25;
        colors(curv>0,:) = .75;

        % Get index for data above the requested thresh (default = 0) and select
        % those data
        if (strcmp(prfParams{p}, 'eccentricity') || strcmp(prfParams{p}, 'eccentricity_smoothed'))
            ii = find(curvPrfData(varExpl.(hemi{h})>0.1));
        end
       
            
        Z = curvPrfData(ii);

        % Convert to 1-266
        Z_ind = round(sz.*((Z-min(Z)) ./ (max(Z)-min(Z))))+1;

        % overlay in colors variable
        colors(ii,:) = cmap(Z_ind,:);

%         colors(~nanIdx) = curvPrfData(~nanIdx);
        
        % Render the triangle mesh
        figure; tH = trimesh(faces, x,y,z);
        
        % Make it look nice
        set(tH, 'LineStyle', 'none', 'FaceColor', 'interp', 'FaceVertexCData',colors)
        axis equal off; colormap(cmap); set(gca, 'CLim', clim)
        colorbar;
        
        % Lighting to make it look glossy
        light('Position',100*[0 1 1],'Style','local')
        
        %     lighting gouraud
        
        % Which mesh are we plotting?
        title(sprintf('Freesurfer mesh: %s %s', hemi{h}, prfParams{p}))
        
        % Print lateral view
        if hemi{h} =='rh'
            vw = [-183.1000   90.0000];
        elseif  hemi{h} =='lh'
            vw = [3.3000  -87.6000];
        end
        set(gca, 'View', vw);
        print(gcf, '-dpng', fullfile(saveFigFolder,   sprintf('%s_FSmesh_%s_%s_lateral', subject, hemi{h}, prfParams{p})))
        
         % Print medial view
        if hemi{h} =='rh'
            vw = [ 0.1000  -90.0000];
        elseif  hemi{h} =='lh'
            vw = [-177.1000   79.6000];
        end
        set(gca, 'View', vw);
        print(gcf, '-dpng', fullfile(saveFigFolder,sprintf('%s_FSmesh_%s_%s_medial', subject, hemi{h}, prfParams{p})))

    end
end

%% 3. Check original versus smoothed parameters (x,y,sigma, beta) on BS surface

thresh = 0;
clims  = [];

for p = 1:length(prfParams)
    
    fprintf('Loading prf data: %s on Brainstorm mesh\n', prfParams{p});
    
    prfDataPthBS = fullfile(prfDataDir, 'surface', 'brainstorm', sprintf('pial.%s', prfParams{p}));
    [curvPrfBS, fnum] = read_curv(prfDataPthBS);
    nanIdx = isnan(curvPrfBS);
    
    cmap = hsv(266);
    
    ttl = sprintf('Brainstorm prf data: %s', prfParams{p});
    visualizeBrainstormMesh(brainstormAnatDir, curvPrfBS, cmap, thresh, clims, [], ttl)
    axis off;
    print(gcf, '-dpng', fullfile(saveFigFolder,sprintf('%s_BSmesh_%s', subject, prfParams{p})))

    
end

%% (2) Visualize Wang atlas on BS surface

fprintf('Loading Wang atlas ROIs on Brainstorm mesh:\n');

wangAtlasPthBS = fullfile(roiDataDir, 'surface','brainstorm', 'pial.all_rois');
[curvWangAtlas, fnum] = read_curv(wangAtlasPthBS);

cmap = hsv(266);

ttl = 'Brainstorm Wang atlas ROIs';
visualizeBrainstormMesh(brainstormAnatDir, curvWangAtlas, cmap, [], [], [], ttl)
colorbar;
print(gcf, '-dpng', fullfile(saveFigFolder,sprintf('%s_BSmesh_wangAtlas', subject)))

