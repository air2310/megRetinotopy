
fs_roi_vz = '/mnt/storage_2/MEG/Retinotopy/Quality_check/wlsubj004/rois/surface/freesurfer/FS_visualize';

surfaces_to_load = {'lh.white','rh.white'};

for n = 1:length(surfaces_to_load)
    
    cur_surf = fullfile(FS_surface,surfaces_to_load{n});
    
    tmp = strsplit(surfaces_to_load{n},'.');
    cur_hs = tmp{1};
    fprintf('Exporting parameters for %s hemisphere:\n',cur_hs);
    
    
    % load our rois which are now as morphometry file formats (saved using
    % write_curv function)
    roi_name = roi_names;
    
    for roi_idx = 1:length(roi_names)
        
        cur_roi = roi_name{roi_idx};
        
        [~,r_name] = fileparts(roi_names{n});
        
        
        
        tmp_roi = read_curv(fullfile(fs_roi_dir,strcat(cur_hs,'.',cur_roi(2:end-4))));
        lindex = find(~isnan(tmp_roi));
        
        %create a label file on the terminal as lh.V1.label for eg:
        labelfile = fullfile(fs_roi_vz,strcat(cur_hs,'.',cur_roi(2:end-4),'.label')); % 'lh.V1.label';
        
        npoints1 = length(lindex);
        npoints = npoints1;
        lxyz = zeros(npoints,3);
        lvals  = zeros(npoints,1);
        
        % open as an ascii file
        fid = fopen(labelfile, 'w') ;
        if(fid == -1)
            fprintf('ERROR: could not open %s\n',labelfile);
            return;
        end
        
        
        % Make sure they are npoints by 1 %
        lindex = reshape(lindex,[npoints 1]);
        lxyz   = reshape(lxyz,[npoints 3]);
        lvals  = reshape(lvals,[npoints 1]);
        
        l = [lindex lxyz lvals];
        fprintf(fid,'%d %f %f %f %f\n',l') ;
        
        fclose(fid) ;
    end
end


%% Saving prf data for visualizing on freesurfer surface

fs_prf_data_vz = '/mnt/storage_2/MEG/Retinotopy/Quality_check/wlsubj004/prf_data/surface/freesurfer/FS_visualize';
   
surfaces_to_load = {'lh.white','rh.white'};

for n = 1:length(surfaces_to_load)
    
    cur_surf = fullfile(FS_surface,surfaces_to_load{n});
    
    tmp = strsplit(surfaces_to_load{n},'.');
    cur_hs = tmp{1};
    fprintf('Exporting parameters for %s hemisphere:\n',cur_hs);
   
    
    %lh_files = dir(fullfile(fs_prf_data,'lh.*'));
    h_files = dir(strcat(fs_prf_data,'/',cur_hs,'.*'));
    
    for nn = 1:length(lh_files)
        
        cur_h_file = h_files(nn).name;
        [~,~,par_name] = fileparts(cur_h_file);
        tmp_pm = cur_h_file; %'lh.polar_angle';
        
        tmp_pm_fpath = fullfile(fs_prf_data,tmp_pm);
        tmp_pm_file = read_curv(tmp_pm_fpath);
        idx_nan = find(isnan(tmp_pm_file));
        tmp_pm_file(idx_nan) = 0;
        
        write_curv(fullfile(fs_prf_data_vz,tmp_pm),tmp_pm_file,0);
        
        fprintf('\n Max of %s = %f \n ',cur_h_file,max(tmp_pm_file));
        fprintf('\n Min of %s = %f \n ',cur_h_file,min(tmp_pm_file));
        
    end
end


%% CHECKING SOME STUFF (delete later)
read_lbl = 0;
if read_lbl == 1
        
    fpath = '/mnt/storage_2/MEG/Retinotopy/Quality_check/wlsubj004/rois/surface/freesurfer/FS_visualize';
    fname = 'lh.V1_check.label';
    roi_file = fullfile(fpath,fname);
    
    % open it as an ascii file
    fid = fopen(roi_file, 'r') ;
    if(fid == -1)
        fprintf('ERROR: could not open %s\n',fname);
        return;
    end
    
    fgets(fid);
    if(fid == -1)
        fprintf('ERROR: could not open %s\n',fname);
        return;
    end
    
    line = fgets(fid) ;
    nv = sscanf(line, '%d') ;
    l = fscanf(fid, '%d %f %f %f %f\n') ;
    l = reshape(l, 5, nv) ;
    l = l' ;
    
    fclose(fid) ;    
   
    
    labelfile = '/mnt/storage_2/MEG/Retinotopy/Quality_check/wlsubj004/rois/surface/freesurfer/FS_visualize/lh.V1_check.label';
      % open as an ascii file
        fid = fopen(labelfile, 'w') ;
        if(fid == -1)
            fprintf('ERROR: could not open %s\n',labelfile);
            return;
        end
        
        lindex = l(:,1);
        npoints1 = length(lindex);
        npoints = npoints1;
        lxyz = zeros(npoints,3);
        lvals  = zeros(npoints,1);
        
        % Make sure they are npoints by 1 %
        lindex = reshape(lindex,[npoints 1]);
        lxyz   = reshape(lxyz,[npoints 3]);
        lvals  = reshape(lvals,[npoints 1]);
        
        l = [lindex lxyz lvals];
        fprintf(fid,'%d %f %f %f %f\n',l') ;
        
        fclose(fid) ;
    
    
    
    
    
end




