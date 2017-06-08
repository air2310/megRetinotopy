function mprf__init_directory_structure
% Not intended to be called by the user
% Function that creates the mprf folder structure.

main_dir = mprf__get_directory('main_dir');

get_dir_for = {
    'rm_stimulus'
    'rm_model'
    'vista_t1'
    'fs_volume'
    'fs_surface'
    'bs_headmodel'
    'bs_anat'
    'fs_surf_roi'
    'bs_surf_roi'
    'fs_surf_prf_data'
    'bs_surf_prf_data'
    'prf_nifti'
    'prf_mat'
    'meg_data'
    'meg_stim'
    'mri_stim'
    'meg_stim_params'
    'lh_rois'
    'rh_rois'
    'meg_imported_stim'
    'model_predictions'
    'model_results'
    };


for n = 1:length(get_dir_for)
    rel_path = mprf__get_directory(get_dir_for{n});
    
    if exist(fullfile(main_dir, rel_path),'dir')        
        fprintf('Directory already exists: %s\n', rel_path);
        
    else
        
        make_dir = fullfile(main_dir, rel_path);
        
        succes = mkdir(make_dir);
        if succes
            fprintf('Created directory %s\n',rel_path)
        end
        
        
    end
    
    mprf__add_relative_path(get_dir_for{n},rel_path, true)
    
end

end

%%%%%
%  if ~isempty(ls(fullfile(main_dir, rel_path)))
%             if delete_all
%                 mprf__empty_dir(fullfile(main_dir, rel_path))
%                 
%             else
%                 
%                 answer = questdlg(sprintf('Directory %s is not empty. Skip or delete?',rel_path),...
%                     'Directory is not empty','Skip','Delete','Delete all','Skip');
%                 
%                 if strcmpi(answer,'skip')
%                     continue
%                 elseif strcmpi(answer,'delete')
%                     mprf__empty_dir(fullfile(main_dir, rel_path));
%                     
%                     
%                 elseif strcmpi(answer,'delete all')
%                     
%                     mprf__empty_dir(fullfile(main_dir, rel_path))
%                     delete_all = true;
%                     
%                 else
%                     return
%                     
%                 end
%             end
%             
%         end
