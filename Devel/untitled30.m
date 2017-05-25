function mprf__export_prf_data_to_fs_surface



%% 
% Do all the necessary preparation, load the mprfSESSION, and check for the
% files we need:
load(fullfile(pwd,'mprfSESSION.mat'));
global mprfSESSION

cd(mprfSESSION.orig.vista_dir);

if exist('mprf_fs_meshFromSurface','file')

else
    addpath(mprfSESSION.root.root_path);
    if exist('mprf_fs_meshFromSurface','file')
        
    else
        error('Could not add root directory to path')
        
    end
end
    


if exist('mrmMapVerticesToGray','file');
    
else
    try
        tbUse('vistasoft');
    catch ME
        error('Vistasoft is not on your Matlab path')
        
    end
    
end




if ~isfield(mprfSESSION,'orig') || ~isfield(mprfSESSION.orig,'vista_dir')
    error('Cannot find mrVista directory path, mprfSession is not properly initialized')
end

if ~isfield(mprfSESSION,'prf_exp') || ~isfield(mprfSESSION.prf_exp,'data_dir') || ...
        ~exist(mprfSESSION.prf_exp.data_dir,'file');
    
    error('No exported pRF parameters have been found. Run mprfSessionSmoothExportPRFParams.m first');
end

if ~isfield(mprfSESSION,'source') || ~isfield(mprfSESSION.source,'fs_surface_dir') || ...
        ~exist(mprfSESSION.source.fs_surface_dir,'dir') || isempty(dir(mprfSESSION.source.fs_surface_dir));
    
    error('No imported freesurfer surfaces found,  mprfSession is not properly initialized');
end





% Open hidden volume
hvol = initHiddenGray;
mmPerVox = viewGet(hvol,'mmpervox');

% Which free surfer surface, stick to the white surface for now, as it
% matches mrVista's layer 1 best (i think)
surfaces_to_load = {'lh.white','rh.white'};


% mprfSessionSmoothExportPRFParams exports the pRF parameters as a
% vector in the same format as they are stored in the Gray view. This way,
% every parameter maps to the corresponding node in the Gray view. We load
% these parameters here:
data_file = dir(fullfile(mprfSESSION.prf_exp.data_dir,'*.mat'));
load(fullfile(mprfSESSION.prf_exp.data_dir,data_file.name));

% Loop over all the parameters stored in the exported data file:
par_names = fieldnames(prf_par_exp);

% Keep track of the X and Y variables, both smoothed and unsmoothed to
% compute the eccenricity and polar angle maps as well:
has_x = false;
has_y = false;
has_x_sm = false;
has_y_sm = false;

% Loop over the surfaces:
for n = 1:length(surfaces_to_load)
    
    cur_surf = fullfile(mprfSESSION.source.fs_surface_dir,surfaces_to_load{n});
    tmp = strsplit(surfaces_to_load{n},'.');
    cur_hs = tmp{1};
    fprintf('Exporting parameters for %s hemisphere:\n',cur_hs)
    
    % Load mesh using fs_meshFromSurface, this creates a mrVista compatible
    % mesh. Using 'my own' function that skips the smoothing:
    mrv_msh = mprf_fs_meshFromSurface(cur_surf);
    fnum = numel(mrv_msh.triangles);
    
    % compute mapping using mrmMapVerticesToGray (mrmMapGrayToVertices):
    cur_v2gmap = mrmMapVerticesToGray(mrv_msh.vertices, viewGet(hvol,'nodes'),...
        mmPerVox);
    
    good_mapping = cur_v2gmap > 0;
    
    % cur_v2gmap is really just an index variable into the gray nodes
    % that tells which gray node is closest to which mesh node. As it uses
    % nearpoints (i.e. nearest neighbour interpolation), I suspect this only
    % works for mrVista based meshes.
    % Would be nice to do something more intelligent than nearest neighbour
    % interpolation, like checking other layers of the Gray graph as well.
    
    % And use the cur_v2gmap to select the pRF parameters of the Gray nodes
    % that are closest to the mesh nodes/vertices. However, cur_v2gmap contains
    % zeros, so we need to deal with them:
    
    for nn = 1:length(par_names)
        
        cur_par_name = par_names{nn};
        tmp = nan(size(cur_v2gmap));
        
        if sum(double(cur_par_name)) == sum(double('beta'))
            tmp_data = squeeze(prf_par_exp.(cur_par_name)(:,:,1));
            tmp(good_mapping) = tmp_data(cur_v2gmap(good_mapping));
            
        else
            tmp(good_mapping) = prf_par_exp.(cur_par_name)(cur_v2gmap(good_mapping));
            
        end
        
        if regexp(cur_par_name,'x\>');
            has_x = true;
            x_pos = tmp;
            disp(cur_par_name)
        elseif regexp(cur_par_name,'y\>');
            has_y = true;
            y_pos = tmp;
            disp(cur_par_name)
        elseif regexp(cur_par_name,'x_smoothed\>');
            has_x_sm = true;
            x_pos_sm = tmp;
            disp(cur_par_name)
        elseif regexp(cur_par_name,'y_smoothed\>');
            has_y_sm = true;
            y_pos_sm = tmp;
            disp(cur_par_name)
            
            
        end
        
        if has_x && has_y
            [tmp_ang, tmp_ecc] = cart2pol(x_pos, y_pos);
            
            fname = mprfExportDataPath('prf_fs_surface',[cur_hs '.polar_angle']);
            write_curv( fname,tmp_ang, fnum);
            fprintf('Polar angle\n')
            
            fname = mprfExportDataPath('prf_fs_surface',[cur_hs '.eccentricity']);
            write_curv( fname,tmp_ecc, fnum);
            fprintf('Eccentricity\n')
            has_x = false;
            has_y = false;
            
        end
        
        
        if has_x_sm && has_y_sm
            [tmp_ang, tmp_ecc] = cart2pol(x_pos_sm, y_pos_sm);
            
            
            fname = mprfExportDataPath('prf_fs_surface',[cur_hs '.polar_angle_smoothed']);
            write_curv( fname,tmp_ang, fnum);
            fprintf('Polar angle smoothed\n')
            
            
            fname = mprfExportDataPath('prf_fs_surface',[cur_hs '.eccentricity_smoothed']);
            write_curv( fname,tmp_ecc, fnum);
            fprintf('Eccentricity smoothed\n')
            
            has_x_sm = false;
            has_y_sm = false;
            
        end
        
        fname = mprfExportDataPath('prf_fs_surface',[cur_hs '.' cur_par_name]);
        write_curv(fname,tmp, fnum);
        fprintf('%s\n', cur_par_name)
        
    end
    
    
    
    
end


cd(mprfSESSION.init.main_dir);
save(fullfile(mprfSESSION.init.main_dir,'mprfSESSION'),'mprfSESSION')
mrvCleanWorkspace;
clear mprfSESSION



