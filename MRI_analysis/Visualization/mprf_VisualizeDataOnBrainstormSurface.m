function bs_msh = mprf_VisualizeDataOnBrainstormSurface(dirPth,cur_surf,saveDir)
% load mesh using fs_meshFromSurface and t_meshFromFreesurfer

% Load the brainstorm surface file:
surf_path = fullfile(dirPth.bs.anatPth);
bs_surf_files = dir(fullfile(surf_path,'*tess_cortex_pial_low.mat'));
bs_msh = mprfMeshFromBrainstorm(fullfile(bs_surf_files.folder,bs_surf_files.name));

% Get prf parameters saved as nifti's
prfParams = {'eccentricity', 'eccentricity_smoothed', 'polar_angle', 'polar_angle_smoothed'...
    'sigma', 'sigma_smoothed', 'varexplained', 'beta', 'recomp_beta', 'wang2015_atlas', 'mask'};


for ii = 1:length(prfParams)
    
    fprintf('(%s):  Visualizing %s on brainstorm surface \n', mfilename, prfParams{ii})
    
    cur_param = strcat(cur_surf,'.',prfParams{ii});
    
    if strcmpi(prfParams{ii},'wang2015_atlas')
        data_in = read_curv(fullfile(dirPth.fmri.saveDataPth_roiBS, cur_param));
    else
        data_in = read_curv(fullfile(dirPth.fmri.saveDataPth_prfBS, cur_param));
    end
    
    if strcmpi(prfParams{ii},'wang2015_atlas') || strcmpi(prfParams{ii},'mask')
        data_in(data_in == 0) = nan;
    end
    
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
    
    % Get list of viewpoints and meshes when saving different images
    viewList={'back','bottom','top','left','right'};
    viewVectors={[pi/2 -pi -pi/2],[pi 0 0],[0 0 pi],[pi/2 -pi/2 0],[-pi/2 -pi/2 0]};
    
    %     viewList={'top'};
    %     viewVectors={[pi/2 -pi -pi/2]};
    
    switch prfParams{ii}
        case {'polar_angle', 'polar_angle_smoothed'}
            diameter = 256;
            mp = [hsv(128); 0.5 0.5 0.5];
            totAngle = 360;
            [img,cmap] = wedgeMap(mp,diameter,totAngle);
            figure(99); imagesc(img); colormap(cmap);
            
        otherwise
            cmap = hsv(256);
    end
    
    % visualize parameter on the mesh
    bs_msh = mprfSessionColorMesh(bs_msh,data_in,cmap,drange,mask);
    for thisView=1:length(viewList)
        cam.actor=0;
        cam.rotation = rotationMatrix3d(viewVectors{thisView});
        mrMesh('localhost',bs_msh.id,'set',cam);
        
        fH = figure('Color', 'w'); clf;
        imagesc(mrmGet(bs_msh, 'screenshot')/255); axis image; axis off;
        print(fH, fullfile(saveDir,sprintf('%s_%s',prfParams{ii},viewList{thisView})), '-dpng');
    end
    
    
    
end

end