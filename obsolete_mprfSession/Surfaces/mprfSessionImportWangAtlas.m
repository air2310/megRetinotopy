function mprfSessionImportWangAtlas

load(fullfile(pwd,'mprfSESSION.mat'))
global mprfSESSION


cd(mprfSESSION.orig.paths.fs_dir)

[atlas_files, atlas_path] = uigetfile('*.mgz',...
    'Please select atlas files you want to import',...
    'MultiSelect','on');

cd(mprfSESSION.init.main_dir)
fprintf('Importing:\n')




for n = 1:length(atlas_files)
    fprintf('%s\n',atlas_files{n});

    surf_vol = MRIread(fullfile(atlas_path, atlas_files{n}));
    [~,tmp] = fileparts(atlas_files{n});
    [~,whs, surf_name] = fileparts(tmp);
    surf_name(surf_name == '.') = '_';
    
    if strfind(surf_name,'wang2015_atlas')
        vol_data = squeeze(surf_vol.vol);
        
        if ~exist('lh_atlas_mask','var') && ...
                strcmpi(whs,'lh')
            lh_atlas_mask = false(size(vol_data));
            lh_atlas_mask = lh_atlas_mask | vol_data > 0;
            
            fname = mprfExportDataPath('wang_roi_fs_surface',[whs '.wang_atlas_mask']);
            write_curv(fname,lh_atlas_mask,1);
            
        elseif ~exist('rh_atlas_mask','var') && ...
                strcmpi(whs,'rh')
            
            rh_atlas_mask = false(size(vol_data));
            rh_atlas_mask = rh_atlas_mask | vol_data > 0;
            
            fname = mprfExportDataPath('wang_roi_fs_surface',[whs '.wang_atlas_mask']);
            write_curv(fname,rh_atlas_mask,1);
            
        end

        fname = mprfExportDataPath('wang_roi_fs_surface',[whs '.wang_atlas']);
        write_curv(fname,squeeze(surf_vol.vol),1);

    elseif strfind(surf_name,'template_areas');
        vol_data = squeeze(surf_vol.vol);
        
        if ~exist('lh_template_mask_roi','var') && ...
                strcmpi(whs,'lh')
            lh_template_mask_roi = false(size(vol_data));
            lh_template_mask_roi = lh_template_mask_roi | vol_data > 0;
            
            fname = mprfExportDataPath('wang_roi_fs_surface',[whs '.wang_template_mask']);
            write_curv(fname,lh_template_mask_roi,1);
            
        elseif ~exist('rh_template_mask_roi','var') && ...
                strcmpi(whs,'rh')
            
            rh_template_mask_roi = false(size(vol_data));
            rh_template_mask_roi = rh_template_mask_roi | vol_data > 0;
            
            fname = mprfExportDataPath('wang_roi_fs_surface',[whs '.wang_template_mask']);
            write_curv(fname,rh_template_mask_roi,1);
            
            
            
        end

        fname = mprfExportDataPath('wang_roi_fs_surface',[whs '.wang_areas']);
        write_curv(fname,squeeze(surf_vol.vol),1);
        
    elseif strfind(surf_name,'v123roi_predict');
        vol_data = abs(squeeze(surf_vol.vol));
        
        if ~exist('lh_predict_mask_roi','var') && ...
                strcmpi(whs,'lh')
            lh_predict_mask_roi = false(size(vol_data));
            lh_predict_mask_roi = lh_predict_mask_roi | vol_data > 0;
            
            fname = mprfExportDataPath('wang_roi_fs_surface',[whs '.wang_predict_mask']);
            write_curv(fname,lh_predict_mask_roi,1);
            
        elseif ~exist('rh_predict_mask_roi','var') && ...
                strcmpi(whs,'rh')
            
            rh_predict_mask_roi = false(size(vol_data));
            rh_predict_mask_roi = rh_predict_mask_roi | vol_data > 0;
            
            fname = mprfExportDataPath('wang_roi_fs_surface',[whs '.wang_predict_mask']);
            write_curv(fname,rh_predict_mask_roi,1);
            
            
            
        end

        fname = mprfExportDataPath('wang_roi_fs_surface',[whs '.wang_v123']);
        write_curv(fname,squeeze(surf_vol.vol),1);
        
        
    elseif any(strfind(surf_name,'template_angle')) || any(strfind(surf_name,'template_eccen'))
       
        vol_data = squeeze(surf_vol.vol);
        
        if ~exist('lh_template_mask','var') && ...
                strcmpi(whs,'lh')
            lh_template_mask = false(size(vol_data));
            lh_template_mask = lh_template_mask | vol_data > 0;
            
            fname = mprfExportDataPath('wang_prf_fs_surface',[whs '.wang_template_mask']);
            write_curv(fname,lh_template_mask,1);
            
        elseif ~exist('rh_template_mask','var') && ...
                strcmpi(whs,'rh')
            
            rh_template_mask = false(size(vol_data));
            rh_template_mask = rh_template_mask | vol_data > 0;
            
            fname = mprfExportDataPath('wang_prf_fs_surface',[whs '.wang_template_mask']);
            write_curv(fname,rh_template_mask,1);
            
            
            
        end
        
        fname = mprfExportDataPath('wang_prf_fs_surface',[whs '.wang' surf_name]);
        write_curv(fname,vol_data,1);

        
    elseif any(strfind(surf_name,'angle_predict')) || any(strfind(surf_name,'eccen_predict'))
        vol_data = squeeze(surf_vol.vol);
        
        if ~exist('lh_predict_mask','var') && ...
                strcmpi(whs,'lh')
            lh_predict_mask = false(size(vol_data));
            lh_predict_mask = lh_predict_mask | vol_data > 0;
            
            fname = mprfExportDataPath('wang_prf_fs_surface',[whs '.wang_predict_mask']);
            write_curv(fname,lh_predict_mask,1);
            
        elseif ~exist('rh_predict_mask','var') && ...
                strcmpi(whs,'rh')
            
            rh_predict_mask = false(size(vol_data));
            rh_predict_mask = rh_predict_mask | vol_data > 0;
            
            fname = mprfExportDataPath('wang_prf_fs_surface',[whs '.wang_predict_mask']);
            write_curv(fname,rh_predict_mask,1);
            
            
            
        end
        
        fname = mprfExportDataPath('wang_prf_fs_surface',[whs '.wang' surf_name]);
        write_curv(fname,vol_data,1);
        
        
        
        
        
    else
        
        error('Unknown surface file type')
        
        
    end
end

if exist('predict_mask','var')

    
    

fprintf('Done. \n')

end





