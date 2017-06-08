function bs_msh = mprfSessionRenderDataOnBrainstormSurface(bs_msh)

load('mprfSESSION.mat')
global mprfSESSION

if ~exist('bs_msh','var') || isempty(bs_msh)
    
    
    load(mprfSESSION.source.bs_head_model);
    [~, surf_name] = fileparts(bs_model.SurfaceFile);
    
    
    bs_msh = mprfMeshFromBrainstorm(mprfSESSION.source.(surf_name));
    
    fname_parts = strsplit(surf_name,'_');
    
    tmp = {'pial','white','mid'};
    
    file_prefix = tmp{ceil(find(~cellfun(@isempty, [cellfun(@(x) strfind(x,tmp{1}), fname_parts, 'UniformOutput',false) ...
        cellfun(@(x) strfind(x,tmp{2}), fname_parts, 'UniformOutput',false) ...
        cellfun(@(x) strfind(x,tmp{3}), fname_parts, 'UniformOutput',false)])) ...
        ./length(fname_parts))};
    
    bs_msh.surf_type   = [file_prefix '.'];

end

cmap_names = {'parula','jet','hsv','hot','colorcube'};

[data_file, data_path] = uigetfile([bs_msh.surf_type '*'],'Select data to render on Brainstorm surface');
cur_cmap = listdlg('ListString',cmap_names,...
    'PromptString','Please select colormap to use',...
    'SelectionMode','single');

surf_data = read_curv(fullfile(data_path, data_file));

eval(['cmap =' cmap_names{cur_cmap} '(256);'])
data_in = surf_data;

answer = questdlg('Please select type of mask to use','Which mask','ROI','Data','Data');
[~, data_hs] = fileparts(data_file);

if strcmpi(answer,'roi')
    try
       mask = read_curv(fullfile(mprfSESSION.roi_exp.brainstorm,[data_hs '.all_rois_mask']));
       mask = ~mask;
        
    catch
        warning('Could not find ROI mask, switching to data mask');
        answer = 'data';
        
    end
    
    
end

if strcmpi(answer,'data')
    
    mask_file = dir(fullfile(data_path,[data_hs '*mask']));
    if length(mask_file) ==1
        mask = read_curv(fullfile(data_path,mask_file(1).name));
        mask = mask == 0;
        
    else
        warning('Could not find data mask, will not use a mask')
        mask = ones(size(surf_data));
        
    end
end

if any([regexp(data_file,'recomp_beta\>') ...
        regexp(data_file,'mresp_smoothed\>') ...
        regexp(data_file,'beta\>') ...
        regexp(data_file,'mresp\>')])
    
    drange = [min(data_in(mask)) prctile(data_in(mask),98)];
    
    
else
    
    drange = [min(data_in) max(data_in)];
    
end


bs_msh = mprfSessionColorMesh(bs_msh,data_in,cmap,drange,mask);


end























