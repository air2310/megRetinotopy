function params = mprf__init_session_gui(root_field, field_names)

global mprfSESSION %#ok<NUSED>
global mprf_fh

if isempty('mprfSESSION')
    error('mprfSESSION is empty')
end

mprf_fh = figure('Units','Normalized');
ctrl_shift_down = [0 .07 0 0];
params.ui.fg_handle = mprf_fh;

for n = 1:length(field_names)
    cur_field = [root_field '.' field_names{n}];
    
    
    params.ui.(field_names{n}).txt = uicontrol(...
        'Style','text',...
        'Units','Normalized',...
        'String',get_field_titles(field_names{n}),...
        'Position',[.01 .93 .25 .03]-ctrl_shift_down*(n-1));
    
    params.ui.(field_names{n}).txt_box = uicontrol(...
        'Style','edit',...
        'Units','Normalized',...
        'String',init_text_field(eval(['mprfSESSION.' cur_field])),...
        'Position',[0.26 .92 .5 .05]-ctrl_shift_down*(n-1),...
        'Callback',{@set_path_from_text_field, root_field, field_names{n}});
    
    params.ui.(field_names{n}).browse = uicontrol(...
        'Style','Pushbutton',...
        'Units','Normalized',...
        'String','Browse',...
        'Position',[.77 .92 .15 .05]-ctrl_shift_down*(n-1),...
        'Callback',{@get_file_and_set_path, root_field, field_names{n}});
    
    params.ui.(field_names{n}).name = field_names{n};
    
end
params.ui.root_field = root_field;

% DONE!
params.ui.done = uicontrol(...
    'Style','Pushbutton',...
    'Units','Normalized',...
    'String','Done',...
    'Position',[.45 .01 .1 .05],...
    'Callback',@exit_ui);

set(mprf_fh,'UserData',params);
uiwait(params.ui.fg_handle);


end



%% Auxiliary functions:
function set_path_from_text_field(obj,~,root_field, field_name)
global mprfSESSION %#ok<NUSED>

dir_name = get(obj,'String'); %#ok<NASGU>

if exist(dir_name,'file') || exist(dir_name,'dir')
    
    eval(['mprfSESSION.' [root_field '.' field_name ] '= dir_name;']);
    eval(['mprfSESSION.updated.' [root_field '.' field_name ] '= true;']);
    
    fprintf('%s.%s set to: %s\n',root_field, field_name, dir_name);
else
    fprintf('ERROR: %s is not a valid file or directory', dir_name)
end

end



function set_str = init_text_field(dir_name)


if iscell(dir_name)
    
    [~, file_names] = cellfun(@fileparts, dir_name,'UniformOutput' , false);
    set_str = '';
    for f = file_names
        set_str = [set_str ' ' f{1}]; %#ok<AGROW>
    end
    
else
    set_str = dir_name;
end

end

function get_file_and_set_path(~,~,root_field, fname)
global mprfSESSION %#ok<NUSED>
global mprf_fh

params = get(mprf_fh,'UserData');


switch lower(fname)
    
    case 'vista_path'        
        params = get_dir_start_dir(params, root_field, fname, 'vista_path',...
            'Please select subject mrVista directory');

    case 'rm_model'
        params = get_file_cdir(params, root_field, fname, {'rm_model','vista_path'},'*.mat',...
            'Please select retinotopic model',false);
        
    case 'bs_anat'
        params = get_dir_start_dir(params, root_field, fname, {'bs_anat','bs_hm'},...
            'Please select brainstorm anatomy directory');
        
    case 'bs_hm'
        params = get_file_cdir(params, root_field, fname, {'bs_anat','bs_hm'},'*.mat',...
            'Please select  Brainstorm headmodel file',false);
        
    case 'fs_dir'
        params = get_dir_start_dir(params, root_field, fname, 'fs_dir',...
            'Please select subject freesurfer directory');
        
    case 'vista_t1'
        params = get_file_cdir(params, root_field, fname, {'vista_t1','vista_class','vista_path'},'*',...
            'Please select mrVista T1 file',false);
        
    case 'vista_class'
        params = get_file_cdir(params, root_field, fname,  {'vista_class','vista_t1','vista_path'},'*',...
            'Please select mrVista T1 class file',false);
        
    case 'lh_rois'
        params = get_file_cdir(params, root_field, fname, {'lh_rois', 'rh_rois','vista_path'},'*.mat',...
            'Please select LEFT HS rois',true);
        
    case 'rh_rois'
        params = get_file_cdir(params, root_field, fname, {'rh_rois', 'lh_rois','vista_path'},'*.mat',...
            'Please select RIGHT HS rois',true);
        
    case 'meg_data'
        params = get_file_cdir(params, root_field, fname, {'meg_data','meg_stim','meg_stim_params'},'*.sqd',...
            'Please select MEG data file',false);
        
    case 'meg_stim'
        params = get_dir_start_dir(params, root_field, fname, {'meg_stim','meg_stim_params','meg_data'},...
            'Please select MEG stimulus directory');
        
    case 'meg_stim_params'
        params = get_dir_start_dir(params, root_field, fname, {'meg_stim_params','meg_stim','meg_data'},...
            'Please select MEG stimulus parameter directory');
        
    case 'mri_stim'
        params = get_dir_start_dir(params, root_field, fname, 'mri_stim',...
            'Please select MEG stimulus directory');
        
    otherwise
        warning('Unknown field name %s',fname)
        
        
end


set(mprf_fh, 'UserData',params)


end

function params = get_dir_start_dir(params, root_field, field_name, check_field,ui_str)

global mprfSESSION %#ok<NUSED>


if ~iscell(check_field)
    check_field = {check_field};
end

for n = 1:length(check_field)
    cur_check_field = check_field{n};

    start_dir = eval(['mprfSESSION.' [root_field '.' cur_check_field] ';']);
    if iscell(start_dir)
        start_dir = start_dir{1};
    end

    if start_dir
        if isdir(start_dir)
            
        else
            start_dir = fileparts(start_dir);
            
        end
        
        break
   
    end
    
end
if start_dir
else
    start_dir = pwd;
end
        

dir_name = uigetdir(start_dir,ui_str);

if dir_name
    params = set_tb_str_and_path(params,root_field, field_name, dir_name);
    
end



end


function params = get_file_cdir(params, root_field, field_name, check_field, file_filt, ui_str,sel_mul_files)

global mprfSESSION %#ok<NUSED>

cur_dir = pwd;

if ~exist('sel_mul_files','var') || isempty(sel_mul_files)
    sel_mul_files = false;
end

if sel_mul_files
    msel = 'On';
elseif ~sel_mul_files
    msel = 'Off';
end

if ~iscell(check_field)
    check_field = {check_field};
end

for n = 1:length(check_field)
    cur_check_field = check_field{n};

    cd_dir = eval(['mprfSESSION.' [root_field '.' cur_check_field] ';']);
    if iscell(cd_dir)
        cd_dir = cd_dir{1};
    end
    
    if cd_dir
        if isdir(cd_dir)
            
        else
            cd_dir = fileparts(cd_dir);
            
        end
        
        break
    end
    
end

if cd_dir
else
    cd_dir = pwd;
end

cd(cd_dir)

[file_name, file_path] = uigetfile(file_filt,ui_str,'MultiSelect',msel);

if ischar(file_name)
    sel_mul_files = false;
end

if (ischar(file_name) && ischar(file_path) && ~sel_mul_files) || ...
    (iscell(file_name) && ischar(file_path) && sel_mul_files)
    params = set_tb_str_and_path(params,root_field, field_name, fullfile(file_path, file_name), sel_mul_files);
    
elseif file_name == 0 && file_path == 0
    
    return
    
else
    error('Unexpected input')
    
    
end

cd(cur_dir)


end

function params = set_tb_str_and_path(params,root_field, field_name, dir_name,mult_files)
global mprfSESSION %#ok<NUSED>

if ~exist('mult_files','var') || isempty(mult_files)
    mult_files = false;
end


if mult_files
    if iscell(dir_name)
        
        [~, file_names] = cellfun(@fileparts, dir_name,'UniformOutput' , false);
        set_str = '';
        for f = file_names
            set_str = [set_str ' ' f{1}]; %#ok<AGROW>
        end
        
        
    else
        error('Expected dir_name to be cell')
    end

else
    if iscell(dir_name)
        error('Dir_name is cell, expected char')
    elseif ischar(dir_name)
        
        set_str = dir_name;
        
    else
        error('Dir_name is not char')
        
    end
    
end
set(params.ui.(field_name).txt_box,...
    'String',set_str);

if iscell(dir_name)
    if all(cellfun(@(x) exist(x, 'file'), dir_name) | cellfun(@(x) exist(x, 'dir'), dir_name))
        eval(['mprfSESSION.' [root_field '.' field_name ] '= dir_name;']);
        eval(['mprfSESSION.updated.' [root_field '.' field_name ] '= true;']);
        fprintf('%s.%s set to: %s\n',root_field, field_name, set_str);
    else
        
        fprintf('ERROR: %s not all inputs are a valid file or directory', set_str)
        
    end
    
elseif ischar(dir_name)
    
    
    if exist(dir_name,'file') || exist(dir_name,'dir')
        eval(['mprfSESSION.' [root_field '.' field_name ] '= dir_name;']);
        eval(['mprfSESSION.updated.' [root_field '.' field_name ] '= true;']);
        fprintf('%s.%s set to: %s\n',root_field, field_name, set_str);
        
    else
        fprintf('ERROR: %s is not a valid file or directory', dir_name)
        
    end
end

end


function exit_ui(~, ~)

params = get(gcf,'UserData');

close(params.ui.fg_handle);

end




function title_name = get_field_titles(field_name)

switch lower(field_name)
    
    case  'vista_path'
        title_name = 'Subject mrVista directory';
        
        
    case 'vista_class'
        title_name = 'Vista class file';
        
    case 'vista_t1'
        title_name = 'Vista anatomy';
        
    case 'bs_hm'
        title_name = 'Brainstorm headmodel';
        
    case 'bs_anat'
        title_name = 'Brainstorm anatomy directory';
        
    case 'fs_dir'
        title_name = 'Subject freesurfer directory';
        
    case 'rm_model'
        title_name = 'Retinotopic model';
        
    case 'lh_rois'
        title_name = 'LEFT hemishpere ROIs';
        
    case 'rh_rois'
        title_name = 'RIGH hemishpere ROIs';
        
    case 'meg_data'
        title_name = 'MEG data';
        
    case 'meg_stim'
        title_name = 'MEG stimulus folder';
        
    case 'meg_stim_params'
        title_name = 'MEG parameter folder';
        
    case 'mri_stim'
        title_name = 'MRI stimulus folder';
        
        
    otherwise
        warning('Unknown field name %s', field_name)
        title_name = '';
        
        
        
        
end







end




