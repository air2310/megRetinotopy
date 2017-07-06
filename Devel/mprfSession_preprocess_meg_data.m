

global mprfSESSION

if isempty('mprfSESSION')
    load('mprfSESSION.mat')
    
end

main_dir = mprf__get_directory('main_dir');
meg_ddir = fullfile(main_dir, mprf__get_directory('meg_data'));
syn_ddir = fullfile(main_dir, mprf__get_direchtory('syn_data'));

has_meg_file = false;
has_syn_file = false;

if exist(meg_ddir,'dir') 
    in_meg_dir = dir(fullfile(meg_ddir, '*.sqd'));
    if ~isempty(in_meg_dir)
        has_meg_file = true;
        
    end
end

if exist(syn_ddir,'dir') 
    in_syn_dir = dir(fullfile(meg_ddir, '*.sqd'));
    if ~isempty(in_syn_dir)
        has_syn_file = true;
        
    end
end

if has_meg_file || has_syn_file

else
    fprintf('Error, no data found')
    return
end










































