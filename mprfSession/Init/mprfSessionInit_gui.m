
function mprfSessionInit_gui
%% GUI Layout:
%
global paths
params.ui.fg_handle = figure('Units','Normalized');

ctrl_shift_down = [0 .08 0 0];


% Defaults:
if ~isfield(paths.orig,'vista_dir')
    paths.orig.vista_dir = '';
end

if ~isfield(paths.orig,'rm_file')
    paths.orig.rm_file = '';
    
end
if ~isfield(paths.orig,'bs_anat_dir')
    paths.orig.bs_anat_dir = '';
    
end
if ~isfield(paths.orig,'bs_model_file')
    paths.orig.bs_model_file = '';
    
end
if ~isfield(paths.orig,'fs_subject_dir')
    paths.orig.fs_subject_dir = '';
    
end

if ~isfield(paths.orig,'mrvt1')
    paths.orig.mrvt1 = '';
    
end


if ~isfield(paths.orig,'mrv_class')
    paths.orig.mrv_class = '';
    
end

% mrVista directory:

params.ui.vista_path.txt = uicontrol(...
    'Style','text',...
    'Units','Normalized',...
    'String','Subject mrVista directory:',...
    'Position',[.01 .93 .25 .03]);

params.ui.vista_path.txt_box = uicontrol(...
    'Style','edit',...
    'Units','Normalized',...
    'String',paths.orig.vista_dir,...
    'Position',[0.26 .92 .5 .05],...
    'Callback',@update_paths);

params.ui.vista_path.browse = uicontrol(...
    'Style','Pushbutton',...
    'Units','Normalized',...
    'String','Browse',...
    'Position',[.77 .92 .15 .05],...
    'Callback',@get_file_and_set_path);



% Retinotopy model file to import:
n=3;
params.ui.rm_path.txt = uicontrol(...
    'Style','text',...
    'Units','Normalized',...
    'String','Retinotopy model file to use:',...
    'Position',[.01 .93 .25 .03]-ctrl_shift_down*n);

params.ui.rm_path.txt_box = uicontrol(...
    'Style','edit',...
    'Units','Normalized',...
    'String',paths.orig.rm_file,...
    'Position',[0.26 .92 .5 .05]-ctrl_shift_down*n,...
    'Callback',@update_paths);

params.ui.rm_path.browse = uicontrol(...
    'Style','Pushbutton',...
    'Units','Normalized',...
    'String','Browse',...
    'Position',[.77 .92 .15 .05]-ctrl_shift_down*n,...
    'Callback',@get_file_and_set_path);


% Brainstorm subject folder:
n=4;
params.ui.bs_path.txt = uicontrol(...
    'Style','text',...
    'Units','Normalized',...
    'String','Subject brainstorm directory:',...
    'Position',[.01 .93 .25 .03]-ctrl_shift_down*n);

params.ui.bs_path.txt_box = uicontrol(...
    'Style','edit',...
    'Units','Normalized',...
    'String',paths.orig.bs_anat_dir ,...
    'Position',[0.26 .92 .5 .05]-ctrl_shift_down*n,...
    'Callback',@update_paths);

params.ui.bs_path.browse = uicontrol(...
    'Style','Pushbutton',...
    'Units','Normalized',...
    'String','Browse',...
    'Position',[.77 .92 .15 .05]-ctrl_shift_down*n,...
    'Callback',@get_file_and_set_path);


% Subject's Brainstorm forward model folder:
n=5;
params.ui.bs_model.txt = uicontrol(...
    'Style','text',...
    'Units','Normalized',...
    'String','Brain storm head model file:',...
    'Position',[.01 .93 .25 .03]-ctrl_shift_down*n);

params.ui.bs_model.txt_box = uicontrol(...
    'Style','edit',...
    'Units','Normalized',...
    'String',paths.orig.bs_model_file,...
    'Position',[0.26 .92 .5 .05]-ctrl_shift_down*n,...
    'Callback',@update_paths);

params.ui.bs_model.browse = uicontrol(...
    'Style','Pushbutton',...
    'Units','Normalized',...
    'String','Browse',...
    'Position',[.77 .92 .15 .05]-ctrl_shift_down*n,...
    'Callback',@get_file_and_set_path);



% Subject's free surfer directory:
n=6;
params.ui.fs_dir.txt = uicontrol(...
    'Style','text',...
    'Units','Normalized',...
    'String','Freesurfer directory:',...
    'Position',[.01 .93 .25 .03]-ctrl_shift_down*n);

params.ui.fs_dir.txt_box = uicontrol(...
    'Style','edit',...
    'Units','Normalized',...
    'String',paths.orig.fs_subject_dir,...
    'Position',[0.26 .92 .5 .05]-ctrl_shift_down*n,...
    'Callback',@update_paths);

params.ui.fs_dir.browse = uicontrol(...
    'Style','Pushbutton',...
    'Units','Normalized',...
    'String','Browse',...
    'Position',[.77 .92 .15 .05]-ctrl_shift_down*n,...
    'Callback',@get_file_and_set_path);

% Subject's mrVista T1:
n=1;
params.ui.mrvt1.txt = uicontrol(...
    'Style','text',...
    'Units','Normalized',...
    'String','mrVista T1:',...
    'Position',[.01 .93 .25 .03]-ctrl_shift_down*n);

params.ui.mrvt1.txt_box = uicontrol(...
    'Style','edit',...
    'Units','Normalized',...
    'String',paths.orig.mrvt1,...
    'Position',[0.26 .92 .5 .05]-ctrl_shift_down*n,...
    'Callback',@update_paths);

params.ui.mrvt1.browse = uicontrol(...
    'Style','Pushbutton',...
    'Units','Normalized',...
    'String','Browse',...
    'Position',[.77 .92 .15 .05]-ctrl_shift_down*n,...
    'Callback',@get_file_and_set_path);

% Subject's mrVista class file:
n=2;
params.ui.mrv_class.txt = uicontrol(...
    'Style','text',...
    'Units','Normalized',...
    'String','mrVista CLASS FILE:',...
    'Position',[.01 .93 .25 .03]-ctrl_shift_down*n);

params.ui.mrv_class.txt_box = uicontrol(...
    'Style','edit',...
    'Units','Normalized',...
    'String',paths.orig.mrv_class,...
    'Position',[0.26 .92 .5 .05]-ctrl_shift_down*n,...
    'Callback',@update_paths);

params.ui.mrv_class.browse = uicontrol(...
    'Style','Pushbutton',...
    'Units','Normalized',...
    'String','Browse',...
    'Position',[.77 .92 .15 .05]-ctrl_shift_down*n,...
    'Callback',@get_file_and_set_path);


% DONE!
params.ui.done = uicontrol(...
    'Style','Pushbutton',...
    'Units','Normalized',...
    'String','Done',...
    'Position',[.45 .01 .1 .05],...
    'Callback',@exit_ui);



set(gcf,'UserData',params);


uiwait(params.ui.fg_handle);

end








%% Auxiliary functions:

function update_paths(obj,junk)
global paths

params = get(gcf,'UserData');
tmp = get(obj,'String');

if obj == params.ui.vista_path.txt_box
    paths.orig.vista_dir = tmp;
    
    
elseif obj == params.ui.rm_path.txt_box
    paths.orig.rm_file = tmp;
    
    
elseif obj == params.ui.bs_path.txt_box
    paths.orig.bs_anat_dir = tmp;
    
    
elseif obj == params.ui.bs_model.txt_box
    paths.orig.bs_model_file = tmp;
    
    
elseif obj == params.ui.fs_dir.txt_box
    paths.orig.fs_subject_dir = tmp;
    
    
elseif obj == params.ui.mrvt1.txt_box
    paths.orig.mrvt1 = tmp;
    
elseif obj == params.ui.mrv_class.txt_box
    paths.orig.mrv_class = tmp;
    
else
    error('Unknown field')
    
    
    
end


set(gcf,'UserData',params);

end



function get_file_and_set_path(obj,junk)
global paths
params = get(gcf,'UserData');

if obj == params.ui.vista_path.browse
    
    dir_name = uigetdir(pwd,'Please select subject mrVista directory');
    
    if dir_name
        set(params.ui.vista_path.txt_box,...
            'String',dir_name);
        
        update_paths(params.ui.vista_path.txt_box)
    end
    
elseif obj == params.ui.rm_path.browse
    cur_dir = pwd;
    
    
    if paths.orig.vista_dir
        cd(paths.orig.vista_dir);
    end
    
    [fname, fpath] = uigetfile('.mat','Please select retinotopy model file');
    
    if ischar(fname) && ischar(fpath)
        
        set(params.ui.rm_path.txt_box,...
            'String',fullfile(fpath, fname));
        
        update_paths(params.ui.rm_path.txt_box)
    end
    
    cd(cur_dir);
    
elseif obj == params.ui.bs_path.browse
    
    start_dir = fileparts(paths.orig.bs_model_file);
    if start_dir
    else
        start_dir = pwd;
    end
    
    dir_name = uigetdir(start_dir,'Please select subject Brainstorm anat directory');
    
    if dir_name
        set(params.ui.bs_path.txt_box,...
            'String',dir_name);
        
        update_paths(params.ui.bs_path.txt_box)
    end
    
    
    
elseif obj == params.ui.bs_model.browse
    cur_dir = pwd;
    
    
    if paths.orig.bs_anat_dir
        cd(paths.orig.bs_anat_dir);
    end
    
    [fname, fpath] = uigetfile('.mat','Please select Brainstorm headmodel file');
    
    if ischar(fname) && ischar(fpath)
        
        set(params.ui.bs_model.txt_box,...
            'String',fullfile(fpath, fname));
        
        update_paths(params.ui.bs_model.txt_box)
    end
    
    cd(cur_dir);
    
elseif obj == params.ui.fs_dir.browse
    
    dir_name = uigetdir(pwd,'Please select subject freesurfer directory');
    
    if dir_name
        set(params.ui.fs_dir.txt_box,...
            'String',dir_name);
        
        update_paths(params.ui.fs_dir.txt_box)
    end
    
    
elseif obj == params.ui.mrvt1.browse
    cur_dir = pwd;
    
    if paths.orig.vista_dir
        cd(paths.orig.vista_dir);
    end
    
    [fname, fpath] = uigetfile('*','Please select mrVista T1 file');
    
    if ischar(fname) && ischar(fpath)
        
        set(params.ui.mrvt1.txt_box,...
            'String',fullfile(fpath, fname));
        
        update_paths(params.ui.mrvt1.txt_box)
    end
    
    cd(cur_dir);
    
elseif obj == params.ui.mrv_class.browse
    
    cur_dir = pwd;
    
    
    if paths.orig.vista_dir
        cd(paths.orig.vista_dir);
    end
    
    [fname, fpath] = uigetfile('*','Please select mrVista CLASS file');
    
    if ischar(fname) && ischar(fpath)
        
        set(params.ui.mrv_class.txt_box,...
            'String',fullfile(fpath, fname));
        
        update_paths(params.ui.mrv_class.txt_box)
    end
    
    cd(cur_dir);
    
    
end





end


function exit_ui(obj, junk)

params = get(gcf,'UserData');

close(params.ui.fg_handle);






end


