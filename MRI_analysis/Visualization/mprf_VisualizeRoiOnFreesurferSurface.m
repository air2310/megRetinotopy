function fs_msh = mprf_VisualizeRoiOnFreesurferSurface(dirPth,cur_surf,saveDir)
% load mesh using fs_meshFromSurface and t_meshFromFreesurfer
% surfaces_to_load = {'lh.pial','lh.white','rh.pial','rh.white'};
    
freesurfer_surface = dirPth.fs.surfPth;
fs_msh = fs_meshFromSurface(fullfile(freesurfer_surface,cur_surf)); 

% And smooth and color it:
fs_msh.smooth_iterations = 128;
fs_msh.smooth_relaxation = 1;
fs_msh = meshSmooth(fs_msh);
fs_msh = meshColor(fs_msh);

% Get prf parameters saved as nifti's
% prfParams = {'eccentricity', 'eccentricity_smoothed', 'polar_angle', 'polar_angle_smoothed'...
%     'sigma', 'sigma_smoothed', 'varexplained', 'beta', 'recomp_beta'};

prfParams = {'V1','V2d','V2v'};

for ii = 1:length(prfParams)
    
    fprintf('(%s):  Visualizing %s on freesurfer surface \n', mfilename, prfParams{ii})
    
    cmap = hsv(256);
    
    cur_hs_tmp = strsplit(cur_surf,'.');
    cur_hs = cur_hs_tmp{1};
    cur_param = strcat(cur_hs,'.',prfParams{ii});
    
    surf_data = read_curv(fullfile(dirPth.fmri.saveDataPth_roiFS, cur_param));
    data_in = surf_data;
    
    
    % mask - usually all the vertices that has a value for the data point.
    mask = true(size(surf_data));
    mask(isnan(surf_data)) = 0;
    
    drange = [min(data_in) max(data_in)];
    
    fs_msh = mprfSessionColorMesh(fs_msh,data_in,cmap,drange,mask);
    
    % Get list of viewpoints and meshes when saving different images
    if strcmpi(cur_hs,'lh')
        viewList={'right'};
        viewVectors={[0 0 pi]};
    elseif strcmpi(cur_hs,'rh')
        viewList={'left'};
        viewVectors={[pi 0 0]};
    end
    
    for thisView=1:length(viewList)
        cam.actor=0;
        cam.rotation = rotationMatrix3d(viewVectors{thisView});
        mrMesh('localhost',fs_msh.id,'set',cam);
        
        fH = figure('Color', 'w'); clf;
        imagesc(mrmGet(fs_msh, 'screenshot')/255); axis image; axis off;
        print(fH, fullfile(saveDir,sprintf('%s_%s',prfParams{ii},viewList{thisView})), '-dpng');
        
    end
    
    

end

