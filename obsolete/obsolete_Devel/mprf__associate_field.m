function field_name = mprf__associate_field(data_type)


switch lower(data_type)

    case 'rm_stimulus'
        field_name = 'source.rm_stim';
        
    case 'rm_model'
        field_name = 'source.rm_model';
        
    case {'vista_t1','vista_class'}
        field_name = 'source.mrv_anat';
        
    case 'fs_volume'
        field_name =  'source.fs_vol';
        
    case 'fs_surface'
        field_name = 'source.fs_surf';
        
    case 'bs_headmodel'
        field_name = 'source.bs_hm';
        
    case 'bs_anat'
        field_name = 'source.bs_anat';
        
    case 'fs_surf_roi'
        field_name = 'rois.surf.fs';
        
    case 'bs_surf_roi'
        field_name = 'rois.surf.bs';
        
        
    case 'fs_surf_prf_data'
        field_name = 'prf_data.surf.fs';
        
        
    case 'bs_surf_prf_data'
        field_name = 'prf_data.surf.bs';
        
        
    case 'prf_nifti'
        field_name = 'prf_data.nifti';
        
    case 'prf_mat'
        field_name = 'prf_data.exp_dir';
        
    case 'meg_data'
        field_name = 'meg.data';
        
    case 'meg_stim'
        field_name = 'meg.stim';
        
    case 'meg_stim_params'
        field_name = 'meg.stim_params';
        
    case 'mri_stim'
        field_name = 'mri.stim';
        
    case 'lh_rois'
        field_name = 'source.vista_rois.lhs';
        
    case 'rh_rois'
        field_name = 'source.vista_rois.rhs';
        
    case 'meg_imported_stim'
        field_name = 'meg.imported_stim';
        
    case 'meg_grid'
        field_name = 'meg.grid';
        
    case 'vista_t1_file'
        field_name = 'source.mrv_t1';
        
    case 'vista_class_file'
        field_name = 'source.mrv_class';
        
    case 'rm_stim_file'
        field_name = 'source.rm_stim_file';
        
    case 'prf_export_mat'
        field_name = 'source.prf_data_file';
        
    case 'bs_headmodel_file'
        field_name = 'source.bs_hm_file';
        
    case 'bs_anat_file'
        field_name = 'source.bs_anat_file';
        
    case 'fs_tag_to_idx'
        field_name = 'rois.surf.fs_tag';
        
    case 'bs_tag_to_idx'
        field_name = 'rois.surf.bs_tag';
        
    case 'model_predictions'
        field_name = 'model.predictions';
        
    case 'model_results'
        field_name = 'model.results';
        
    otherwise
        error('Unkown data type')
        
        
        
        
        
        
        
        
        
        
        
        
        
        
end