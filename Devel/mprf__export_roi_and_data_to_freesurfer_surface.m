function success = mprf__export_roi_and_data_to_freesurfer_surface(w_exp)

% check for mprfSESSION
global mprfSESSION
if isempty(mprfSESSION)
    error('mprfSESSION is empty')
end

% Get main directory
main_dir = mprf__get_directory('main_dir');

success = 0;

% Check if we have everything we need to do what we want:
if strcmpi(w_exp, 'both')
    % We need freesurfer surfaces (both LH and RH)
    % We need exported prf parameters (.mat file)
    % We need the ROIs, separately for left and right
    % We need a working mrVista session directory:
    if mprfSESSION.has.fs_surface.lh_white && ...
            mprfSESSION.has.fs_surface.rh_white && ...
            mprfSESSION.has.prf_exported && ...
            mprfSESSION.has.lh_rois && mprfSESSION.has.rh_rois && ...
            mprfSESSION.has.vista_path
        
    else
        fprintf(['Error: do not have the required data to export pRF parameters'...
            'to freesurfer surfaces.\n Please update session and import the'...
            'required data and export pRF parameters'])
        return
        
        
        
        
    end
    
    
elseif strcmpi(w_exp, 'prf_data')
    if mprfSESSION.has.fs_surface.lh_white && ...
            mprfSESSION.has.fs_surface.rh_white && ...
            mprfSESSION.has.prf_exported && ...
            mprfSESSION.has.vista_path
        
    else
        fprintf(['Error: do not have the required data to export pRF parameters'...
            'to freesurfer surfaces.\n Please update session and import the'...
            'required data and export pRF parameters'])
        return
        
        
        
        
    end
    
    
elseif strcmpi(w_exp, 'rois')
    if mprfSESSION.has.fs_surface.lh_white && ...
            mprfSESSION.has.fs_surface.rh_white && ...
            mprfSESSION.has.lh_rois && mprfSESSION.has.rh_rois && ...
            mprfSESSION.has.vista_path
        
    else
        fprintf(['Error: do not have the required data to export pRF parameters'...
            'to freesurfer surfaces.\n Please update session and import the'...
            'required data and export pRF parameters'])
        return
        
        
        
        
    end
    
end

% Navigate to mrVista session:
cd(mprfSESSION.orig.paths.vista_path);

% Check for neccessary files:
if exist('mprf_fs_meshFromSurface','file')
    
else
    addpath(mprfSESSION.root.root_path);
    if exist('mprf_fs_meshFromSurface','file')
        
    else
        error('Could not add root directory to path')
        
    end
end
% Again, better to make this one call to tbUse
if exist('mrmMapVerticesToGray','file');
    
else
    try
        tbUse({'retMEG','vistasoft','knk'});
    catch ME
        error('Vistasoft is not on your Matlab path')
        
    end
    
end


% Open hidden volume
hvol = initHiddenGray;
mmPerVox = viewGet(hvol,'mmpervox');


success_rois = false;
success_data = false;
% Do all the work:
if strcmpi(w_exp, 'both') || strcmpi(w_exp, 'prf_data')
    success_data = mprf__export_prf_data_to_fs_surface(hvol, mmPerVox, {'lh.white','rh.white'});
    
end


if strcmpi(w_exp, 'both') || strcmpi(w_exp, 'rois')
    success_rois = mprf__export_vista_rois_to_fs_surface(hvol, mmPerVox, {'lh.white','rh.white'});
    
    
end

if success_rois && success_data
    success = 3;
    
elseif success_rois && ~success_data
    success = 1;
    
elseif success_data && ~success_rois
    success = 2;
    
end


cd(main_dir)
mrvCleanWorkspace;




end

















