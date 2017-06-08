function [yes_to_all, succes] = mprf__import_data(orig_path, file_type, out_file_name, yes_to_all)

if ~exist('yes_to_all','var') || isempty(yes_to_all)
    yes_to_all = false;
end

overwrite = false;

if ischar(orig_path)
    dest_path = fullfile(mprf__get_directory('main_dir'), mprf__get_directory(file_type),out_file_name);
    
    sys_cmnd = sprintf('cp %s %s',orig_path,dest_path);
    
    if exist(dest_path,'file')
        if yes_to_all
            overwrite = true;
            
        else
            answer = questdlg(sprintf('%s is already imported, overwrite?',out_file_name),...
                'Overwrite','No','Yes','Yes to all','No');
            
            switch lower(answer)
                
                case 'no'
                    
                    
                case 'yes'
                    overwrite = true;
                    
                case 'yes to all'
                    overwrite = true;
                    yes_to_all = true;
                    
            end
            
        end
        
        if overwrite
            if system(sys_cmnd);
                tmp = strsplit(orig_path, ' ');
                orig_path_02 = tmp{1};
                
                
                for n = 2:length(tmp)
                    orig_path_02 = [orig_path_02 '\ ' tmp{n}];
                end
                
                sys_cmnd = sprintf('cp %s %s',orig_path_02,dest_path);
                
                
                if system(sys_cmnd)
                    
                    warning('Could not copy %s \n', file_type);
                else
                    fprintf('Imported %s \n', file_type)
                end
                
                
            else
                fprintf('Imported %s \n', file_type)
                
            end
            
        end
        
        
    else
        if system(sys_cmnd);
            tmp = strsplit(orig_path, ' ');
            orig_path_02 = tmp{1};
            
            
            for n = 2:length(tmp)
                orig_path_02 = [orig_path_02 '\ ' tmp{n}];
            end
            
            sys_cmnd = sprintf('cp %s %s',orig_path_02,dest_path);
            
            
            if system(sys_cmnd)
                
                warning('Could not copy %s \n', file_type);
            else
                fprintf('Imported %s \n', file_type)
            end
            
            
        else
            fprintf('Imported %s \n', file_type)
            
        end
    end
    
else
    warning('No valid file name or path obtained, skipping')
end

succes = exist(dest_path,'file');






end


