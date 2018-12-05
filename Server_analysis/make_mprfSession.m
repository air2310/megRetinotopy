function mprfSESSION = make_mprfSession(sub_sess_dir)

mprfSESSION.init.main_dir = sub_sess_dir;

% source directories
mprfSESSION.source.rm_stim = strcat(sub_sess_dir,'/source/mrvista/rm_stimulus');
mprfSESSION.source.rm_model = strcat(sub_sess_dir,'/source/mrvista/ret_model');
mprfSESSION.source.mrv_anat = strcat(sub_sess_dir,'/source/mrvista/anatomy');
mprfSESSION.source.fs_vol = strcat(sub_sess_dir,'/source/freesurfer/volumes');
mprfSESSION.source.fs_surf = strcat(sub_sess_dir,'/source/freesurfer/surfaces');
mprfSESSION.source.bs_hm = strcat(sub_sess_dir,'/source/brainstorm/head_model');
mprfSESSION.source.bs_anat = strcat(sub_sess_dir,'/source/brainstorm/anatomy');
mprfSESSION.source.vista_rois.lhs = strcat(sub_sess_dir,'/source/mrvista/rois/lhs');
mprfSESSION.source.vista_rois.rhs = strcat(sub_sess_dir,'/source/mrvista/rois/rhs');
mprfSESSION.source.mrv_t1 = strcat(sub_sess_dir,'/source/mrvista/anatomy/t1.nii.gz');
mprfSESSION.source.mrv_class = strcat(sub_sess_dir,'/source/mrvista/anatomy/t1_class.nii.gz');
mprfSESSION.source.bs_anat_file = strcat(sub_sess_dir,'/source/brainstorm/anatomy/subjectimage_T1.mat');

d = dir(fullfile(sub_sess_dir,'/source/mrvista/rm_stimulus/*.mat'));
mprfSESSION.source.rm_stim_file = fullfile(d.folder, d.name);
mprfSESSION.source.prf_data_file = strcat(sub_sess_dir,'/prf_data/data_dir/exported_prf_params.mat');

% Surface ROI directories
mprfSESSION.rois.surf.fs = strcat(sub_sess_dir,'/rois/surface/freesurfer');
mprfSESSION.rois.surf.bs = strcat(sub_sess_dir,'/rois/surface/brainstorm');
mprfSESSION.rois.surf.fs_tag = strcat(sub_sess_dir,'/rois/surface/freesurfer/fs_tag_to_idx');
mprfSESSION.rois.surf.bs_tag = strcat(sub_sess_dir,'/rois/surface/brainstorm/bs_tag_to_idx');

% prf data
mprfSESSION.prf_data.surf.fs = strcat(sub_sess_dir,'/prf_data/surface/freesurfer');
mprfSESSION.prf_data.surf.bs = strcat(sub_sess_dir,'/prf_data/surface/brainstorm');
mprfSESSION.prf_data.nifti = strcat(sub_sess_dir,'/prf_data/nifti');
mprfSESSION.prf_data.exp_dir = strcat(sub_sess_dir, '/prf_data/data_dir');

% MEG data
mprfSESSION.meg.data = strcat(sub_sess_dir,'/data/meg/raw');
mprfSESSION.meg.stim = strcat(sub_sess_dir,'/stimuli/meg/stimulus_files');
mprfSESSION.meg.stim_params = strcat(sub_sess_dir, '/stimuli/meg/param_files');
mprfSESSION.meg.imported_stim = strcat(sub_sess_dir,'/stimuli/meg/imported_stimulus');

% MRI stimulus
mprfSESSION.mri = strcat(sub_sess_dir,'/stimuli/mri/stimulus_files');

%
mprfSESSION.model.predictions =  strcat(sub_sess_dir,'/modeling/predictions');
mprfSESSION.model.results =  strcat(sub_sess_dir,'/modeling/results');

save(fullfile(sub_sess_dir,'mprfSESSION.mat'),'mprfSESSION');

end