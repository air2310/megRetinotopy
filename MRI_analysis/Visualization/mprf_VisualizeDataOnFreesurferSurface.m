function fs_msh = mprf_VisualizeDataOnFreesurferSurface(dirPth,cur_surf,saveDir)
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
prfParams = {'eccentricity', 'eccentricity_smoothed', 'polar_angle', 'polar_angle_smoothed'...
  'sigma', 'sigma_smoothed', 'varexplained', 'beta', 'recomp_beta'};


for ii = 1:length(prfParams)
    
    fprintf('(%s):  Visualizing %s on freesurfer surface \n', mfilename, prfParams{ii})
   
    cur_hs_tmp = strsplit(cur_surf,'.');
    cur_hs = cur_hs_tmp{1};
    cur_param = strcat(cur_hs,'.',prfParams{ii});
   
    switch prfParams{ii}
        case {'polar_angle', 'polar_angle_smoothed'}
            if strcmpi(cur_hs,'lh')
                load(which('WedgeMapLeft_pRF.mat'));
                cmap = modeInformation.cmap;
                cmap = cmap(65:end,:);
            elseif strcmpi(cur_hs,'rh')
                load(which('WedgeMapRight_pRF.mat'))
                cmap = modeInformation.cmap;
                cmap = cmap(65:end,:);
            end
        otherwise
            cmap = hsv(256);
    end
    
    
    surf_data = read_curv(fullfile(dirPth.fmri.saveDataPth_prfFS, cur_param));
    data_in = surf_data;
    
    
    % mask - usually all the vertices that has a value for the data point.
    mask = true(size(data_in));
    mask(isnan(data_in)) = 0;
    
    if any([regexp(prfParams{ii},'recomp_beta\>') ...
            regexp(prfParams{ii},'mresp_smoothed\>') ...
            regexp(prfParams{ii},'beta\>') ...
            regexp(prfParams{ii},'mresp\>')])
        
        data_in(data_in == 0) = nan;
        drange = [min(data_in(mask)) prctile(data_in(mask),95)];
        
        
    else
        
        drange = [min(data_in) max(data_in)];
        
    end
    
    
    
    fs_msh = mprfSessionColorMesh(fs_msh,data_in,cmap,drange,mask);
    
    % Get list of viewpoints and meshes when saving different images
    % viewList={'back','left','right','bottom','top'};
    % viewVectors={[pi -pi/2 0],[pi 0 0],[0 0 pi],[pi/2 -pi/2 0],[-pi/2 -pi/2 0]};
    
    viewList={'back','left'};%
    viewVectors={[pi -pi/2 0],[pi 0 0]};
    
    for thisView=1:length(viewList)
        cam.actor=0;
        cam.rotation = rotationMatrix3d(viewVectors{thisView});
        mrMesh('localhost',fs_msh.id,'set',cam);
        
        fH = figure('Color', 'w'); clf;
        imagesc(mrmGet(fs_msh, 'screenshot')/255); axis image; axis off;
        %print(fH, fullfile(saveDir,sprintf('%s_%s',prfParams{ii},viewList{thisView})), '-dpng');
        saveas(fH, fullfile(saveDir,sprintf('%s_%s',prfParams{ii},viewList{thisView})), 'png');
    end
    
    

end

end