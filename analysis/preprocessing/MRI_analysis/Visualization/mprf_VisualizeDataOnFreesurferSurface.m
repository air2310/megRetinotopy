function fs_msh = mprf_VisualizeDataOnFreesurferSurface(dirPth,cur_surf,saveDir, opt)
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
  'sigma', 'sigma_smoothed', 'beta', 'recomp_beta'};

if opt.mri.useHCPAveMaps
    dirPth.fmri.saveDataPth_prfFS = fullfile(dirPth.fmri.saveDataPth,'HCP7TAveRetinotopyPred','prfFS');   
    prfParams = {'eccentricity', 'polar_angle', 'sigma',  'beta'};
end

fH = figure('Color', 'w'); clf;
for ii = 1:length(prfParams)
    
    fprintf('(%s):  Visualizing %s on freesurfer surface \n', mfilename, prfParams{ii})
   
    cur_hs_tmp = strsplit(cur_surf,'.');
    cur_hs = cur_hs_tmp{1};
    
    ve = read_curv(fullfile(dirPth.fmri.saveDataPth_prfFS, [cur_hs '.varexplained']));
    vemask = ve>opt.mri.varExplThresh(1);
    eccen = read_curv(fullfile(dirPth.fmri.saveDataPth_prfFS, [cur_hs '.eccentricity']));
    eccenmask = eccen<opt.mri.eccThresh(2);
    
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
    if opt.mri.useHCPAveMaps
        cmap = hsv(256);
    end
    
    surf_data = read_curv(fullfile(dirPth.fmri.saveDataPth_prfFS, cur_param));
    surf_data(~vemask) = NaN;
    surf_data(~eccenmask) = NaN;
    data_in = surf_data;
    
    % mask - usually all the vertices that has a value for the data point.
    mask = true(size(data_in));
    mask(isnan(data_in)) = false;

    cAxisLim = [min(data_in) max(data_in)];
        
    fs_msh = mprfSessionColorMesh(fs_msh,data_in,cmap,cAxisLim,mask);
    
    % Get list of viewpoints and meshes when saving different images
    viewList={'back','left','right','bottom','top'};
    viewVectors={[pi -pi/2 0],[pi 0 0],[0 0 pi],[pi/2 -pi/2 0],[-pi/2 -pi/2 0]};
    
    for thisView=1:length(viewList)
        cam.actor=0;
        cam.rotation = rotationMatrix3d(viewVectors{thisView});
        mrMesh('localhost',fs_msh.id,'set',cam);
        
        figure(fH); cla;
        imagesc(mrmGet(fs_msh, 'screenshot')/255); axis image; axis off;
        cmapbar = cmap; caxis(cAxisLim); colormap(cmapbar); colorbar; 
        
        title(sprintf('%s %s %s',cur_hs, prfParams{ii},viewList{thisView})); drawnow;
        saveas(fH, fullfile(saveDir,sprintf('%s_%s_%s', cur_hs, prfParams{ii},viewList{thisView})), 'png');
    end
    
    

end

end