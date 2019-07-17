function bs_msh = mprf_VisualizeDataOnBrainstormSurface(dirPth,cur_surf,saveDir)
% load mesh using fs_meshFromSurface and t_meshFromFreesurfer
% surfaces_to_load = {'lh.pial','lh.white','rh.pial','rh.white'};

% Load the brainstorm surface file:
surf_path = fullfile(dirPth.bs.anatPth);
bs_surf_files = dir(fullfile(surf_path,'*tess_cortex_pial_low.mat'));
bs_msh = mprfMeshFromBrainstorm(fullfile(bs_surf_files.folder,bs_surf_files.name));

%cmap_names = {'parula','jet','hsv','hot','colorcube','phase_left_prf','phase_right_prf'};
%cur_cmap = listdlg('ListString',cmap_names,'PromptString','Please select colormap to use','SelectionMode','single');

% Get prf parameters saved as nifti's
prfParams = {'eccentricity', 'eccentricity_smoothed', 'polar_angle', 'polar_angle_smoothed'...
   'sigma', 'sigma_smoothed', 'varexplained', 'beta', 'recomp_beta','wang2015_atlas','mask'};

%prfParams = {'polar_angle_smoothed'};

hemispheres = {'lh','rh'};

for ii = 1:length(prfParams)
    
    fprintf('(%s):  Visualizing %s on brainstorm surface \n', mfilename, prfParams{ii})   
    
    cur_param = strcat(cur_surf,'.',prfParams{ii});
    
    surf_data = read_curv(fullfile(dirPth.fmri.saveDataPth_prfBS, cur_param));
    data_in = surf_data;
    
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
%    viewList={'back','left','right','bottom','top'};
%    viewVectors={[pi -pi/2 0],[pi 0 0],[0 0 pi],[pi/2 -pi/2 0],[-pi/2 -pi/2 0]};
    
    viewList={'top'};
    viewVectors={[pi/2 -pi -pi/2]};
    
    switch prfParams{ii}
        case {'polar_angle', 'polar_angle_smoothed'}
            load(which('WedgeMapLeft_pRF.mat'));
            cmap_lh = modeInformation.cmap;
            cmap_lh = cmap_lh(65:end,:);
            load(which('WedgeMapRight_pRF.mat'))
            cmap_rh = modeInformation.cmap;
            cmap_rh = cmap_rh(65:end,:);
            
             
            for cur_hs = 1:length(hemispheres)
                % visualize parameter on the mesh
                if strcmpi(hemispheres{cur_hs},'lh')
                    cmap = cmap_lh;
                elseif strcmpi(hemispheres{cur_hs},'rh')
                    cmap = cmap_rh;
                end
                
                bs_msh = mprfSessionColorMesh(bs_msh,data_in,cmap,drange,mask);
                for thisView=1:length(viewList)
                    cam.actor=0;
                    cam.rotation = rotationMatrix3d(viewVectors{thisView});
                    mrMesh('localhost',bs_msh.id,'set',cam);
                    
                    fH = figure('Color', 'w'); clf;
                    imagesc(mrmGet(bs_msh, 'screenshot')/255); axis image; axis off;
                    print(fH, fullfile(saveDir,sprintf('%s_%s_%s',prfParams{ii},viewList{thisView},hemispheres{cur_hs})), '-dpng');
                end
            end
            
        otherwise
            cmap = hsv(256);
            
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

end