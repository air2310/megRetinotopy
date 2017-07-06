function varargout = mprf__get_directory(varargin)
% Idea is to get the directory for the data type in varargin

global mprfSESSION
varargout = cell(size(varargin));


for n = 1:length(varargin)
    
    switch lower(varargin{n})
        
        
        case 'rm_stimulus'
            varargout{n} = 'source/mrvista/rm_stimulus';
            
        case 'rm_model'
            varargout{n} = 'source/mrvista/ret_model';
            
        case {'vista_t1','vista_class'}
            varargout{n} = 'source/mrvista/anatomy';
            
        case 'fs_volume'
            varargout{n} =  'source/freesurfer/volumes';
            
        case 'fs_surface'
            varargout{n} = 'source/freesurfer/surfaces';
            
        case {'bs_headmodel','bs_headmodel_file'}
            varargout{n} = 'source/brainstorm/head_model';
            
        case {'bs_anat','bs_anat_file'}
            varargout{n} = 'source/brainstorm/anatomy';
            
        case 'fs_surf_roi'
            varargout{n} = 'rois/surface/freesurfer';
            
        case 'bs_surf_roi'
            varargout{n} = 'rois/surface/brainstorm';
            
            
        case 'fs_surf_prf_data'
            varargout{n} = 'prf_data/surface/freesurfer';
            
            
        case 'bs_surf_prf_data'
            varargout{n} = 'prf_data/surface/brainstorm';
            
        case 'prf_nifti'
            varargout{n} = 'prf_data/nifti';
            
        case 'prf_mat'
            varargout{n} = 'prf_data/data_dir';
            
        case 'main_dir'
            varargout{n} = mprfSESSION.init.main_dir;
            
        case 'meg_data'
            varargout{n} = 'data/meg/raw';
            
        case {'meg_stim','meg_grid'}
            varargout{n} = 'stimuli/meg/stimulus_files';
            
        case 'meg_stim_params'
            varargout{n} = 'stimuli/meg/param_files';
            
        case 'mri_stim'
            varargout{n} = 'stimuli/mri/stimulus_files';
            
        case 'lh_rois'
            varargout{n} = 'source/mrvista/rois/lhs';
            
        case 'rh_rois'
            varargout{n} = 'source/mrvista/rois/rhs';
            
        case 'meg_imported_stim'
            varargout{n} = 'stimuli/meg/imported_stimulus';
            
        case 'model_predictions'
            varargout{n}  = 'modeling/predictions';
            
        case 'model_results'
            varargout{n}  = 'modeling/results';

        case 'syn_data'
            varargout{n}  = 'modeling/synthetic';
            
        case 'meg_preproc'
            varargout{n}  = 'data/meg/preproc';
            
        case 'syn_preproc'
            varargout{n}  = 'modeling/synthetic/preproc';
            
        otherwise
            
            fprintf('ERROR: unknow type %s\n',varargin{n})
            print_list = {
                'meg_preproc'
                'syn_data'
                'model_results'
                'model_predictions'
                'meg_imported_stim'
                'rh_rois'
                'lh_rois'
                'mri_stim'
                'meg_stim_params'
                'meg_stim'
                'meg_grid'
                'meg_data'
                'main_dir'
                'prf_mat'
                'prf_nifti'
                'bs_surf_prf_data'
                'fs_surf_prf_data'
                'bs_surf_roi'
                'fs_surf_roi'
                'bs_anat'
                'bs_anat_file'
                'bs_headmodel'
                'bs_headmodel_file'
                'fs_surface'
                'fs_volume'
                'vista_t1'
                'vista_class'
                'rm_model'
                'rm_stimulus'};
            
            
            fprintf('Please use any of these keys:\n')
            for nn = length(print_list):-1:1
               fprintf('%s\n',print_list{nn})
                
                
            end
            
            
            
            
            
    end
    
end


end