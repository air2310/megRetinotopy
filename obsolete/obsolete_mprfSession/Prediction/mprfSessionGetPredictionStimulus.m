function stim_file_to_use = mprfSessionGetPredictionStimulus

global mprfSESSION

make_stimulus = true;

if ~isfield(mprfSESSION,'pred') || ~isfield(mprfSESSION.pred,'stimuli')
    
else
    if iscell(mprfSESSION.pred.stimuli) && ~isempty(mprfSESSION.pred.stimuli{1})
        fdir = fileparts(mprfSESSION.pred.stimuli{1});
        stim_files = dir(fullfile(fdir,'*.mat'));
        
        if length(stim_files) == length(mprfSESSION.pred.stimuli)
            % Okay, found as many stimulus files as there are paths
            
        else
            
            
            % Field found, but the amount of paths and stimuli in the
            % directory are not equal
            
            mprfSESSION.pred.stimuli = cell(1,length(stim_files));
            
            for n = 1:length(stim_files)
                mprfSESSION.pred.stimuli{1,n} = fullfile(fdir,stim_files(n).name);
            end
            
            
        end

        
        stim_file_names = {stim_files.name};
        stim_file_names = [stim_file_names 'Make new stimulus'];
        
        answer = listdlg('ListString',stim_file_names,...
            'PromptString','Select existing stimulus, or make new one');
        
        stim_opt = stim_file_names{answer};
        
        if strcmpi(stim_opt,'make new stimulus')
            
        else
            stim_file_to_use = fullfile(fdir,stim_file_names{answer});
            make_stimulus = false;
        end
        
        
    else
        
        
    end
    
    
    
end


if make_stimulus
    answer = questdlg(sprintf('Do you want to use the retinotopic stimulus or an other stimulus?'),...
        'Select stimulus',...
        'Use retinotopic stimulus',...
        'Other stimulus',...
        'Other stimulus');
    
    if strcmpi(answer, 'use retinotopic stimulus')
        
        load(mprfSESSION.source.rm_stim);
        stim_file_to_use = mprfSessionPrepareStimulus(rm_stim);
        
    elseif strcmpi(answer,'other stimulus')
        [s_name, s_path] = uigetfile('*.mat','Please select stimulus file to use');
        cur_dir = pwd;
        cd(s_path);
        [g_name, g_path] = uigetfile('*.mat','Please select corresponding grid file');
        cd(cur_dir);
        
        stim_file = fullfile(s_path, s_name);
        grid_file = fullfile(g_path, g_name);
        
        stim_file_to_use = mprfSessionPrepareStimulus({stim_file, grid_file});
        
    end
    
end










end