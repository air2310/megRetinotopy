function mprf__empty_dir(cur_dir)


indir = dir(cur_dir);

for nn = 1:length(indir)
    if strcmpi(indir(nn).name, '.') || strcmpi(indir(nn).name,'..')
        
    else
        delete(fullfile(cur_dir,indir(nn).name));
        
        if ~exist(fullfile(cur_dir,indir(nn).name),'file')
           fprintf('Deleted %s\n', fullfile(cur_dir,indir(nn).name))
            
        end
        
    end
    
end



end