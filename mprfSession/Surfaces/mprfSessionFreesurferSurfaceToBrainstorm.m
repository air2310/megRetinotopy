function mprfSessionFreesurferSurfaceToBrainstorm(w_export)

% TO DO: when exportin ROIs be sure to check for ROIs that are defined in
% only one of the two hemispheres.

load(fullfile(pwd,'mprfSESSION.mat'));
global mprfSESSION

if strcmpi(w_export,'prf') || strcmpi(w_export,'roi') || strmpci('wang')
else
    error('Unrecognized export, use ''prf'' ,''wang'' or ''roi'' only')
end

if exist('mrmMapVerticesToGray','file');
    
else
    try
        tbUse('vistasoft');
    catch ME
        error('Vistasoft is not on your Matlab path')
        
    end
    
end

if exist('mprf_fs_meshFromSurface','file')
    
else
    addpath(mprfSESSION.root.root_path);
    if exist('mprf_fs_meshFromSurface','file')
        
    else
        error('Could not add root directory to path')
        
    end
end


if ~isfield(mprfSESSION,'source') || ~isfield(mprfSESSION.source,'bs_head_model')  ||  ...
        ~exist(mprfSESSION.source.bs_head_model,'file')
    error('Cannot find Brain storm head model, mprfSession is not properly initialized')
end

if strcmpi(w_export,'prf')
    if ~isfield(mprfSESSION,'prf_exp') || ~isfield(mprfSESSION.prf_exp,'fs_surface_data')  ||  ...
            isempty(dir(mprfSESSION.prf_exp.fs_surface_data))
        error('Cannot find PRF data exported to FreeSurfer surface run mprfSessionExportPRFDataToFreeSurferSurface first')
    end
    
elseif strcmpi(w_export,'roi')
    
    if ~isfield(mprfSESSION,'roi_exp') || ~isfield(mprfSESSION.roi_exp,'freesurfer')  ||  ...
            isempty(dir(mprfSESSION.roi_exp.freesurfer))
        error('Cannot find ROIs exported to FreeSurfer surface run mprfSessionExportVistaROIToFreesurferSurface first')
    end
    
end


if ~isfield(mprfSESSION,'source') || ~isfield(mprfSESSION.source,'fs_surface_dir') || ...
        ~exist(mprfSESSION.source.fs_surface_dir,'dir') || isempty(dir(mprfSESSION.source.fs_surface_dir));
    
    error('No imported freesurfer surfaces found,  mprfSession is not properly initialized');
end

if ~isfield(mprfSESSION,'source') || ~isfield(mprfSESSION.source,'subjectimage_T1') || ...
        ~exist(mprfSESSION.source.subjectimage_T1,'file')
    error('No imported Brainstorm anatomy found,  mprfSession is not properly initialized');
    
end


load(mprfSESSION.source.bs_head_model);
bs_head_model_surf = bs_model.SurfaceFile;
bs_mri = load(mprfSESSION.source.subjectimage_T1);

clear bs_model

hs_to_load = {'lh','rh'};        % In this order as BS concatenates the hemispheres like this


if any(strcmpi('pial',strsplit(bs_head_model_surf,'_')))
    surfaces_to_load = {'pial'};
    
elseif any(strcmpi('white',strsplit(bs_head_model_surf,'_')))
    surfaces_to_load = {'white'};
    
else
    error('Could not determine the surface to load');
    
end


fprintf('Exporting surface data:\n')
for n = 1:length(surfaces_to_load)
    verts = [];
    
    
    for nn = 1:length(hs_to_load)
        cur_surf = [hs_to_load{nn} '.' surfaces_to_load{n}];
        surf_file = fullfile(mprfSESSION.source.fs_surface_dir,cur_surf);
        tmp_verts = mne_read_surface(surf_file); % Use the same routine as Brainstorm uses, otherwise the vertices are slightly off and intersectCols does not find all the matching vertices.
        
        tmp_verts = bsxfun(@plus, tmp_verts, [128 129 128] / 1000);
        tmp_verts = [tmp_verts'; ones(1,size(tmp_verts,1))]; %% Actually verts)
        tmp_verts = [bs_mri.SCS.R, bs_mri.SCS.T./1000; 0 0 0 1] * tmp_verts;
        tmp_verts = tmp_verts(1:3,:)';
        
        verts = [verts; tmp_verts]; %#ok<AGROW>
        
    end
    
    
    eval(['bs_surf_file = mprfSESSION.source.tess_cortex_' surfaces_to_load{n} '_low;'])
    load(bs_surf_file,'Vertices');
    
    [~, bs_vert_idx, bs_vert_idx2] = intersectCols(verts', Vertices');
    
    if numel(bs_vert_idx) == size(Vertices,1);
        
    else
        error('Could only find %d of %d brainstorm vertices in transformed FreeSurfer vertices',...
            length(bs_vert_idx), size(Vertices,1));
    end
    
    if strcmpi(w_export,'prf')
        pname = mprfExportDataPath('prf_fs_surface','');
        
    elseif strcmpi(w_export,'roi')
        pname = mprfExportDataPath('roi_fs_surface','');
        
    end
    
    lh_files = dir(fullfile(pname,'lh.*'));
    
    for nn = 1:length(lh_files)
        
        cur_lh_file = lh_files(nn).name;
        [~,~,par_name] = fileparts(cur_lh_file);
        if strcmpi(par_name,'.mgz')
            
            
        else
            par_name = par_name(2:end);
            cur_rh_file = ['rh.' par_name];
            
            
            both_data = [read_curv(fullfile(pname,cur_lh_file));...
                read_curv(fullfile(pname,cur_rh_file))];
            
            both_bs_data_out = nan(size(bs_vert_idx2));
            both_bs_data_out(bs_vert_idx2) = both_data(bs_vert_idx);
            
            cur_out_file = [surfaces_to_load{n} '.' par_name];
            
            if strcmpi(w_export,'prf')
                fname = mprfExportDataPath('prf_bs_surface',cur_out_file);
                
            elseif strcmpi(w_export,'roi')
                fname = mprfExportDataPath('roi_bs_surface',cur_out_file);
                
            end
            
            write_curv(fname,both_bs_data_out,1);
            fprintf('%s\n',cur_out_file);
        end
        
    end
    
    
    
end

fprintf('Done.\n')


cd(mprfSESSION.init.main_dir);
save(fullfile(mprfSESSION.init.main_dir,'mprfSESSION'),'mprfSESSION')
mrvCleanWorkspace;
clear mprfSESSION




end











