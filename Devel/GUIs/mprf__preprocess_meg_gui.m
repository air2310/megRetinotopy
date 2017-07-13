function varargout = mprf__preprocess_meg_gui(varargin)
% MPRF__PREPROCESS_MEG_GUI MATLAB code for mprf__preprocess_meg_gui.fig
%      MPRF__PREPROCESS_MEG_GUI, by itself, creates a new MPRF__PREPROCESS_MEG_GUI or raises the existing
%      singleton*.
%
%      H = MPRF__PREPROCESS_MEG_GUI returns the handle to a new MPRF__PREPROCESS_MEG_GUI or the handle to
%      the existing singleton*.
%
%      MPRF__PREPROCESS_MEG_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MPRF__PREPROCESS_MEG_GUI.M with the given input arguments.
%
%      MPRF__PREPROCESS_MEG_GUI('Property','Value',...) creates a new MPRF__PREPROCESS_MEG_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before mprf__preprocess_meg_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to mprf__preprocess_meg_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help mprf__preprocess_meg_gui

% Last Modified by GUIDE v2.5 11-Jul-2017 15:58:33

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @mprf__preprocess_meg_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @mprf__preprocess_meg_gui_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before mprf__preprocess_meg_gui is made visible.
function mprf__preprocess_meg_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to mprf__preprocess_meg_gui (see VARARGIN)

% Choose default command line output for mprf__preprocess_meg_gui
handles.output = hObject;
handles.files.meg = varargin{1};
handles.files.stim = varargin{3};
handles.do_syn = varargin{4};

handles.dir.main = mprf__get_directory('main_dir');
handles.dir.meg_data = mprf__get_directory('meg_data');
handles.dir.stim = mprf__get_directory('meg_imported_stim');

if handles.do_syn
   [fpath, fname, fext] = fileparts(varargin{2}{1});
   handles.files.syn = {[fname fext]};
   handles.dir.syn_data = fpath;
   set(handles.uipanel1,'SelectedObject',handles.rb_syn);
   set(handles.lb_load_data,'Value',1);

   
   
else
    handles.files.syn = varargin{2};
    handles.dir.syn_data = mprf__get_directory('syn_data');
    
end


tmp = load(fullfile(handles.dir.main, handles.dir.stim, handles.files.stim{1}));
fname = fieldnames(tmp);
handles.stim = tmp.(fname{1});
handles.stim.full_im = handles.stim.full_im ./ max(handles.stim.full_im(:));
im_range = [1 size(handles.stim.full_im,3)];

set(handles.lbl_stim_max,'String',num2str(im_range(2)));
set(handles.stim_slider,...
    'Min',1,...
    'Max',im_range(2),...
    'SliderStep',[(1/(im_range(2) - im_range(1))) (10/(im_range(2) - im_range(1)))],...
    'Value',1)

axes(handles.axes1)
colormap('gray');

handles = updateStim(hObject, handles);
uipanel1_SelectionChangeFcn(hObject, eventdata, handles);
set(handles.cb_save_at_every_step,'Value',1);
set(handles.cb_hp_filt,'Value',1);
set(handles.cb_pre_proc_data,'Value',1);
set(handles.cb_denoise_data,'Value',1);
set(handles.epoch_lower,'String',num2str(150));
set(handles.epoch_upper,'String',num2str(1100));

set(handles.txt_data_channels,'String','0:156');
set(handles.txt_trigger_channels,'String',' 160:167');
set(handles.txt_diode_channel,'String','191');


% Update handles structure
guidata(hObject, handles);
 
% UIWAIT makes mprf__preprocess_meg_gui wait for user response (see UIRESUME)
uiwait(handles.figure1);


function handles = updateStim(hObject, handles)

cur_val = ceil(get(handles.stim_slider,'Value'));

axes(handles.axes1);
imagesc(handles.stim.full_im(:,:,cur_val), [0 1]);
set(handles.axes1,'YTick',[]);
set(handles.axes1,'XTick',[]);
set(handles.lbl_cur_im, 'String',num2str(cur_val));

% --- Outputs from this function are returned to the command line.
function varargout = mprf__preprocess_meg_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
preproc.do.save = logical(get(handles.cb_save_at_every_step,'Value'));
preproc.do.filt = logical(get(handles.cb_hp_filt ,'Value'));
preproc.do.preproc = logical(get(handles.cb_pre_proc_data ,'Value'));
preproc.do.denoise = logical(get(handles.cb_denoise_data ,'Value'));
preproc.params.epoch.low = str2double(get(handles.epoch_lower,'String'));
preproc.params.epoch.include = str2double(get(handles.epoch_upper,'String'));
tmp = get(handles.lb_load_data,'String');

if get(handles.uipanel1,'SelectedObject') == handles.rb_meg;
    preproc.data.type = 'meg_data';
elseif get(handles.uipanel1,'SelectedObject') == handles.rb_syn;
    preproc.data.type = 'syn_data';
end

preproc.data.file = tmp{get(handles.lb_load_data,'Value')};
preproc.stimulus = handles.stim;

handles.output.channels.data = eval(get(handles.txt_data_channels,'String'));
handles.output.channels.triggers = eval(get(handles.txt_trigger_channels,'String'));
handles.output.channels.diode = eval(get(handles.txt_diode_channel,'String'));


handles.output.preproc = preproc;
% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on slider movement.
function stim_slider_Callback(hObject, eventdata, handles)
% hObject    handle to stim_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

handles = updateStim(hObject, handles);
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function stim_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to stim_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes when selected object is changed in uipanel1.
function uipanel1_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel1 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

if get(handles.uipanel1,'SelectedObject') == handles.rb_meg
    if isempty(handles.files.meg)
        set(handles.lb_load_data,'String','No data found');
    else
        set(handles.lb_load_data,'String',handles.files.meg);
    end
    
elseif get(handles.uipanel1,'SelectedObject') == handles.rb_syn
    if isempty(handles.files.syn)
        set(handles.lb_load_data,'String','No data found');
    else
        set(handles.lb_load_data,'String',handles.files.syn);
    end
end


% --- Executes on selection change in lb_load_data.
function lb_load_data_Callback(hObject, eventdata, handles)
% hObject    handle to lb_load_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns lb_load_data contents as cell array
%        contents{get(hObject,'Value')} returns selected item from lb_load_data


% --- Executes during object creation, after setting all properties.
function lb_load_data_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lb_load_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in cb_save_at_every_step.
function cb_save_at_every_step_Callback(hObject, eventdata, handles)
% hObject    handle to cb_save_at_every_step (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cb_save_at_every_step


% --- Executes on button press in cb_hp_filt.
function cb_hp_filt_Callback(hObject, eventdata, handles)
% hObject    handle to cb_hp_filt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cb_hp_filt


% --- Executes on button press in cb_pre_proc_data.
function cb_pre_proc_data_Callback(hObject, eventdata, handles)
% hObject    handle to cb_pre_proc_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cb_pre_proc_data


% --- Executes on button press in cb_denoise_data.
function cb_denoise_data_Callback(hObject, eventdata, handles)
% hObject    handle to cb_denoise_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cb_denoise_data



function epoch_upper_Callback(hObject, eventdata, handles)
% hObject    handle to epoch_upper (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of epoch_upper as text
%        str2double(get(hObject,'String')) returns contents of epoch_upper as a double


% --- Executes during object creation, after setting all properties.
function epoch_upper_CreateFcn(hObject, eventdata, handles)
% hObject    handle to epoch_upper (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function epoch_lower_Callback(hObject, eventdata, handles)
% hObject    handle to epoch_lower (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of epoch_lower as text
%        str2double(get(hObject,'String')) returns contents of epoch_lower as a double


% --- Executes during object creation, after setting all properties.
function epoch_lower_CreateFcn(hObject, eventdata, handles)
% hObject    handle to epoch_lower (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pb_go.
function pb_go_Callback(hObject, eventdata, handles)
% hObject    handle to pb_go (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

figure1_CloseRequestFcn(hObject, eventdata,handles,true);


% --- Executes on button press in pb_cancel.
function pb_cancel_Callback(hObject, eventdata, handles)
% hObject    handle to pb_cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
figure1_CloseRequestFcn(hObject, eventdata,handles,false);


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles,do_go)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~exist('do_go','var') || isempty(do_go)
    do_go = false;
end

handles.output.do_go = do_go;

guidata(hObject, handles);

% Hint: delete(hObject) closes the figure
if isequal(get(handles.figure1, 'waitstatus'), 'waiting')
    uiresume(handles.figure1);
    
else
    % Hint: delete(hObject) closes the figure
    delete(handles.figure1);
end



function txt_data_channels_Callback(hObject, eventdata, handles)
% hObject    handle to txt_data_channels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_data_channels as text
%        str2double(get(hObject,'String')) returns contents of txt_data_channels as a double


% --- Executes during object creation, after setting all properties.
function txt_data_channels_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_data_channels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txt_trigger_channels_Callback(hObject, eventdata, handles)
% hObject    handle to txt_trigger_channels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_trigger_channels as text
%        str2double(get(hObject,'String')) returns contents of txt_trigger_channels as a double


% --- Executes during object creation, after setting all properties.
function txt_trigger_channels_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_trigger_channels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txt_diode_channel_Callback(hObject, eventdata, handles)
% hObject    handle to txt_diode_channel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_diode_channel as text
%        str2double(get(hObject,'String')) returns contents of txt_diode_channel as a double


% --- Executes during object creation, after setting all properties.
function txt_diode_channel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_diode_channel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
