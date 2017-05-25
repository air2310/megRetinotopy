function success = mprf__export_prf_data_to_fs_surface(hvol,mmPerVox, surfaces_to_load)
success = false;


global mprfSESSION
if isempty(mprfSESSION)
    error('mprfSESSION is empty')
end

if ~exist('surfaces_to_load','var') || isempty(surfaces_to_load)
    surfaces_to_load = {'lh.white','rh.white'};
    
end
main_dir = mprf__get_directory('main_dir');

%%
% Do all the necessary preparation, load the mprfSESSION, and check for the
% files we need:

% Which free surfer surface, stick to the white surface for now, as it
% matches mrVista's layer 1 best (i think)


% mprfSessionSmoothExportPRFParams exports the pRF parameters as a
% vector in the same format as they are stored in the Gray view. This way,
% every parameter maps to the corresponding node in the Gray view. We load
% these parameters here:
data_file = dir(fullfile(main_dir, mprf__get_directory('prf_mat'),'*.mat'));
load(fullfile(main_dir, mprf__get_directory('prf_mat'),data_file.name));

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
    
    cur_surf = fullfile(main_dir, mprf__get_directory('fs_surface'),surfaces_to_load{n});
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
            
            fname = fullfile(main_dir, mprf__get_directory('fs_surf_prf_data'),[cur_hs '.polar_angle']);
            write_curv( fname,tmp_ang, fnum);
            fprintf('Polar angle\n')
            
            fname = fullfile(main_dir, mprf__get_directory('fs_surf_prf_data'),[cur_hs '.eccentricity']);
            write_curv( fname,tmp_ecc, fnum);
            fprintf('Eccentricity\n')
            has_x = false;
            has_y = false;
            
        end
        
        
        if has_x_sm && has_y_sm
            [tmp_ang, tmp_ecc] = cart2pol(x_pos_sm, y_pos_sm);
            
            
            fname = fullfile(main_dir, mprf__get_directory('fs_surf_prf_data'),[cur_hs '.polar_angle_smoothed']);
            write_curv( fname,tmp_ang, fnum);
            fprintf('Polar angle smoothed\n')
            
            
            fname = fullfile(main_dir, mprf__get_directory('fs_surf_prf_data'),[cur_hs '.eccentricity_smoothed']);
            write_curv( fname,tmp_ecc, fnum);
            fprintf('Eccentricity smoothed\n')
            
            has_x_sm = false;
            has_y_sm = false;
            
        end
        
        fname = fullfile(main_dir, mprf__get_directory('fs_surf_prf_data'),[cur_hs '.' cur_par_name]);
        write_curv(fname,tmp, fnum);
        fprintf('%s\n', cur_par_name)
        
    end
        
end

if exist(fname,'file')
    mprfSESSION.has.fs_surf_data = true;
    
else
    mprfSESSION.has.fs_surf_data = false;
    
end
success = true;
cd(main_dir);
save(fullfile(mprfSESSION.init.main_dir,'mprfSESSION'),'mprfSESSION')
mrvCleanWorkspace;

end

