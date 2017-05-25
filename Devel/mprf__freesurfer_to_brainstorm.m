function mprf__freesurfer_to_brainstorm(w_export, bs_model_file)

% TO DO: when exportin ROIs be sure to check for ROIs that are defined in
% only one of the two hemispheres.

global mprfSESSION
if isempty('mprfSESSION')
    error('mprfSESSION is empty')
end

main_dir = mprf__get_directory('main_dir');

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



load(bs_model_file);
bs_head_model_surf = bs_model.SurfaceFile;
bs_mri = load(fullfile(main_dir, mprf__get_file('bs_anat')));

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
        surf_file = fullfile(main_dir,mprf__get_directory('fs_surface'),cur_surf);
        tmp_verts = mne_read_surface(surf_file); % Use the same routine as Brainstorm uses, otherwise the vertices are slightly off and intersectCols does not find all the matching vertices.
        
        tmp_verts = bsxfun(@plus, tmp_verts, [128 129 128] / 1000);
        tmp_verts = [tmp_verts'; ones(1,size(tmp_verts,1))]; %% Actually verts)
        tmp_verts = [bs_mri.SCS.R, bs_mri.SCS.T./1000; 0 0 0 1] * tmp_verts;
        tmp_verts = tmp_verts(1:3,:)';
        
        verts = [verts; tmp_verts]; %#ok<AGROW>
        
    end
    
    surf_path = fullfile(main_dir, mprf__get_directory('bs_anat'));
    
    bs_surf_files = dir(fullfile(surf_path,'*tess_cortex*.mat'));
    this_bs_file = find(~cellfun(@isempty,strfind({bs_surf_files.name},surfaces_to_load{n})));
    
    if isempty(this_bs_file)
        error('Could not find Brainstorm surface file for %s',surfaces_to_load{n})
    end
    load(fullfile(surf_path,bs_surf_files(this_bs_file).name),'Vertices');
    
    [~, bs_vert_idx, bs_vert_idx2] = intersectCols(verts', Vertices');
    
    if numel(bs_vert_idx) == size(Vertices,1);
        
    else
        error('Could only find %d of %d brainstorm vertices in transformed FreeSurfer vertices',...
            length(bs_vert_idx), size(Vertices,1));
    end
    
    if strcmpi(w_export,'prf')
        pname = mprf__get_directory('fs_surf_prf_data');

    elseif strcmpi(w_export,'roi')
        pname = mprf__get_directory('fs_surf_roi');

        
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
                fname = fullfile(mprf__get_directory('bs_surf_prf_data'),cur_out_file);
                
            elseif strcmpi(w_export,'roi')
                fname = fullfile(mprf__get_directory('bs_surf_roi'),cur_out_file);
                
            end
            
            write_curv(fname,both_bs_data_out,1);
            fprintf('%s\n',cur_out_file);
        end
        
    end
    
    
    
end

fprintf('Done.\n')

if strcmpi(w_export,'roi')
    check_dir = 'bs_surf_prf_data';
    
    
elseif strcmpi(w_export,'prf')
    check_dir = 'bs_surf_roi';
    
end


if ~isempty(dir(fullfile(main_dir, mprf__get_directory(check_dir))))
    mprfSESSION.has.(check_dir) = true;
    
else
    mprfSESSION.has.(check_dir) = false;
    
end

cd(main_dir);
save(fullfile(main_dir,'mprfSESSION'),'mprfSESSION')
mrvCleanWorkspace;
clear mprfSESSION



































end





















