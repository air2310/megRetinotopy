function fs_msh = mprfSessionRenderDataOnFreesurferSurface(fs_msh)

load('mprfSESSION.mat')
global mprfSESSION

if ~exist('fs_msh','var') || isempty(fs_msh)
    % Which free surfer surface
    
    cd(mprfSESSION.source.fs_surface_dir)
    
    [surf_name, surf_path] = uigetfile('*','Please select surface');
    
    
    cd(mprfSESSION.init.main_dir)
    
    [~,file_filter] = fileparts(surf_name);
    file_filter = [file_filter ,'.*'];
    
    % load mesh using fs_meshFromSurface and t_meshFromFreesurfer
    % surfaces_to_load = {'lh.pial','lh.white','rh.pial','rh.white'};
    fs_msh = fs_meshFromSurface(fullfile(surf_path,surf_name));
    
    % And smooth and color it:
    fs_msh.smooth_iterations = 128;
    fs_msh.smooth_relaxation = 1;
    fs_msh = meshSmooth(fs_msh);
    fs_msh = meshColor(fs_msh);
    fs_msh.mprf.file_filter = file_filter;
    
end



cmap_names = {'parula','jet','hsv','hot','colorcube','phase_left_prf','phase_right_prf'};
[data_file, data_path] = uigetfile(fs_msh.mprf.file_filter,'Select data to render');
cur_cmap = listdlg('ListString',cmap_names,...
    'PromptString','Please select colormap to use',...
    'SelectionMode','single');

surf_data = read_curv(fullfile(data_path, data_file));


if strcmpi(cmap_names{cur_cmap},'phase_left_prf')
    try
        load(which('WedgeMapLeft_pRF.mat'))
        cmap = modeInformation.cmap;
        
        
    catch
        cmap =hsv(256);
        
        
    end
    
    
elseif strcmpi(cmap_names{cur_cmap},'phase_right_prf')
    try
        load(which('WedgeMapRight_pRF.mat'))
        cmap = modeInformation.cmap;
        
    catch
        cmap =hsv(256);
        
    end
    
else
    eval(['cmap =' cmap_names{cur_cmap} '(256);'])
    
end


data_in = surf_data;
[~, data_hs,data_fname] = fileparts(data_file);


if any(strfind(data_fname,'predict')) || any(strfind(data_fname,'v123'))
    mask_file = fullfile(data_path,[data_hs '.wang_predict_mask']);
    inv_mask = false;
    
elseif any(strfind(data_fname,'template')) || any(strfind(data_fname,'areas'))
    mask_file = fullfile(data_path,[data_hs '.wang_template_mask']);
    inv_mask = false;
    
    
elseif strfind(data_fname,'atlas')
    mask_file = fullfile(data_path,[data_hs '.wang_atlas_mask']);
    inv_mask = false;
    
    
else
    mask_file = fullfile(data_path,[data_hs '.mask']);
    inv_mask = true;
    
    
end

if exist(mask_file,'file')
    
    try
        mask = read_curv(mask_file);
        if inv_mask
            mask = mask == 0;
            
        else 
            mask = logical(mask);
            
        end
    catch
        
        
        warning('Could not load mask file, no masking');
        mask = true(size(surf_data));
        
    end
    
    
else
    mask = true(size(surf_data));
    
end

answer = questdlg('Limite pRF data to ROI extend?');

if strcmpi(answer,'yes')
    [roi_mask_name, roi_mask_path] = uigetfile('*','Please select ROI mask to use');
    
    try
        mask = logical((read_curv(fullfile(roi_mask_path, roi_mask_name)) > 0) .* mask);
        mask = ~mask;
    catch
        fprintf('Could not load ROI mask');
    end
end


if any([regexp(data_file,'recomp_beta\>') ...
        regexp(data_file,'mresp_smoothed\>') ...
        regexp(data_file,'beta\>') ...
        regexp(data_file,'mresp\>')])
    
    data_in(data_in == 0) = nan;
    drange = [min(data_in(mask)) prctile(data_in(mask),95)];
    
    
else
    
    drange = [min(data_in) max(data_in)];
    
end

fs_msh = mprfSessionColorMesh(fs_msh,data_in,cmap,drange,mask);

end









