function out = mprfUpdateSession_Orig_paths

global paths

orig_fields = fieldnames(paths.orig);
out = [];

for n = 1:length(orig_fields)
    
    cur_field = orig_fields{n};
    
    
    if strcmpi(cur_field,'vista_dir')
         if mprfUpdateSession(3, 'add_orig_paths','vista_dir',paths.orig.(cur_field)) == 7;
             out = 1;
         else
             error('Directory %s for the vista directory does not exist or is not a directory',...
                 paths.orig.(cur_field))
         end
        
    elseif strcmpi(cur_field,'fs_subject_dir')
         if mprfUpdateSession(3, 'add_orig_paths','fs_subject_dir',paths.orig.(cur_field)) == 7;
            out = 1;
         else
             error('Directory %s for the freesurfer directory does not exist or is not a directory',...
                 paths.orig.(cur_field))
         end
        
    elseif strcmpi(cur_field, 'rm_file')
        if mprfUpdateSession(3, 'add_orig_paths','ret_model',paths.orig.(cur_field)) == 2;
            out = 1;
        else
            error('File %s for the retinotopy model does not exist or is not a file',...
                paths.orig.(cur_field))
            
        end
             
    elseif strcmpi(cur_field,'bs_anat_dir')
         if mprfUpdateSession(3, 'add_orig_paths','bs_anat_dir',paths.orig.(cur_field)) == 7;
            out = 1;
         else
             warning('Directory %s for the Brainstorm anatomy directory does not exist or is not a directory',...
                 paths.orig.(cur_field))
         end
        
        
    elseif strcmpi(cur_field,'bs_model_file')
         if mprfUpdateSession(3, 'add_orig_paths','bs_head_model',paths.orig.(cur_field)) == 2;
            out = 1;
        else
            warning('File %s for the Bs head model does not exist or is not a file',...
                paths.orig.(cur_field))
            
        end
         
        
    elseif strcmpi(cur_field,'mrvt1')
         if mprfUpdateSession(3, 'add_orig_paths','mr_vista_t1',paths.orig.(cur_field)) == 2;
            out = 1;
        else
            error('File %s for the mrVista T1 does not exist or is not a file',...
                paths.orig.(cur_field))
            
        end
         
        
    elseif strcmpi(cur_field,'mrv_class')
         if mprfUpdateSession(3, 'add_orig_paths','mr_vista_class',paths.orig.(cur_field)) == 2;
            out = 1;
        else
            error('File %s for the mrVista CLASS file does not exist or is not a file',...
                paths.orig.(cur_field))
            
         end
         
    elseif strcmpi(cur_field,'root_path')
        if mprfUpdateSession(3, 'add_root_path','root_path',paths.orig.(cur_field)) == 7;
            out = 1;
        else
            error('Directory %s for the root directory does not exist or is not a directory',...
                paths.orig.(cur_field))
            
         end
        
        
    else
        warning('Unknown field %s, skipping',cur_field);
        
        
        
        
    end
    
end

end








