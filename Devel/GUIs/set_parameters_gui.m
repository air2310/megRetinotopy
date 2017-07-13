function varargout = set_parameters_gui(varargin)
% SET_PARAMETERS_GUI MATLAB code for set_parameters_gui.fig
%      SET_PARAMETERS_GUI, by itself, creates a new SET_PARAMETERS_GUI or raises the existing
%      singleton*.
%
%      H = SET_PARAMETERS_GUI returns the handle to a new SET_PARAMETERS_GUI or the handle to
%      the existing singleton*.
%
%      SET_PARAMETERS_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SET_PARAMETERS_GUI.M with the given input arguments.
%
%      SET_PARAMETERS_GUI('Property','Value',...) creates a new SET_PARAMETERS_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before set_parameters_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to set_parameters_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help set_parameters_gui

% Last Modified by GUIDE v2.5 13-Jul-2017 12:38:22

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @set_parameters_gui_OpeningFcn, ...
    'gui_OutputFcn',  @set_parameters_gui_OutputFcn, ...
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


% --- Executes just before set_parameters_gui is made visible.
function set_parameters_gui_OpeningFcn(hObject, eventdata, handles, varargin) %#ok<*INUSL>
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to set_parameters_gui (see VARARGIN)
handles.ch_model.type = varargin{1}{1};
handles.ch_model.def_params = mprf__get_default_model_parameters(handles.ch_model.type);
set__widget_default_values(hObject, handles, handles.ch_model.def_params);
% Choose default command line output for set_parameters_gui
handles.output = hObject;
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes set_parameters_gui wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = set_parameters_gui_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
out.params  =  handles.ch_model.params_out;
out.type    =  handles.ch_model.type;

if strcmp(get(get(handles.bg_radio_beta,'SelectedObject'),'String'), 'Equal maximum prediction')
    out.params.beta_equal_pred = true;
    out.params.beta_equal_beta = false;
    
elseif strcmp(get(get(handles.bg_radio_beta,'SelectedObject'),'String'), 'Equal beta weight')
    out.params.beta_equal_pred = false;
    out.params.beta_equal_beta = true;
    
end

varargout = {out};
% Get default command line output from handles structure
delete(handles.figure1)


% --- Executes on button press in pb_click_me.
function pb_click_me_Callback(hObject, eventdata, handles) %#ok<*DEFNU>
% hObject    handle to pb_click_me (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
guidata(handles.figure1, handles);


% --- Executes on selection change in pu_prf_x.
function pu_prf_x_Callback(hObject, eventdata, handles)
% hObject    handle to pu_prf_x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pu_prf_x contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pu_prf_x


% --- Executes during object creation, after setting all properties.
function pu_prf_x_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pu_prf_x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in pu_prf_y.
function pu_prf_y_Callback(hObject, eventdata, handles)
% hObject    handle to pu_prf_y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pu_prf_y contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pu_prf_y


% --- Executes during object creation, after setting all properties.
function pu_prf_y_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pu_prf_y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in pu_prf_sigma.
function pu_prf_sigma_Callback(hObject, eventdata, handles)
% hObject    handle to pu_prf_sigma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pu_prf_sigma contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pu_prf_sigma


% --- Executes during object creation, after setting all properties.
function pu_prf_sigma_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pu_prf_sigma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in pu_prf_beta.
function pu_prf_beta_Callback(hObject, eventdata, handles)
% hObject    handle to pu_prf_beta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pu_prf_beta contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pu_prf_beta


% --- Executes during object creation, after setting all properties.
function pu_prf_beta_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pu_prf_beta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function txt_fix_x_Callback(hObject, eventdata, handles)
% hObject    handle to txt_fix_x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_fix_x as text
%        str2double(get(hObject,'String')) returns contents of txt_fix_x as a double


% --- Executes during object creation, after setting all properties.
function txt_fix_x_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_fix_x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txt_x_range_Callback(hObject, eventdata, handles)
% hObject    handle to txt_x_range (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_x_range as text
%        str2double(get(hObject,'String')) returns contents of txt_x_range as a double


% --- Executes during object creation, after setting all properties.
function txt_x_range_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_x_range (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txt_fix_y_Callback(hObject, eventdata, handles)
% hObject    handle to txt_fix_y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_fix_y as text
%        str2double(get(hObject,'String')) returns contents of txt_fix_y as a double


% --- Executes during object creation, after setting all properties.
function txt_fix_y_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_fix_y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txt_y_range_Callback(hObject, eventdata, handles)
% hObject    handle to txt_y_range (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_y_range as text
%        str2double(get(hObject,'String')) returns contents of txt_y_range as a double


% --- Executes during object creation, after setting all properties.
function txt_y_range_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_y_range (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txt_sigma_fix_Callback(hObject, eventdata, handles)
% hObject    handle to txt_sigma_fix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_sigma_fix as text
%        str2double(get(hObject,'String')) returns contents of txt_sigma_fix as a double


% --- Executes during object creation, after setting all properties.
function txt_sigma_fix_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_sigma_fix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txt_sigma_range_Callback(hObject, eventdata, handles)
% hObject    handle to txt_sigma_range (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_sigma_range as text
%        str2double(get(hObject,'String')) returns contents of txt_sigma_range as a double


% --- Executes during object creation, after setting all properties.
function txt_sigma_range_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_sigma_range (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txt_beta_fix_Callback(hObject, eventdata, handles)
% hObject    handle to txt_beta_fix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_beta_fix as text
%        str2double(get(hObject,'String')) returns contents of txt_beta_fix as a double


% --- Executes during object creation, after setting all properties.
function txt_beta_fix_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_beta_fix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txt_beta_range_Callback(hObject, eventdata, handles)
% hObject    handle to txt_beta_range (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_beta_range as text
%        str2double(get(hObject,'String')) returns contents of txt_beta_range as a double


% --- Executes during object creation, after setting all properties.
function txt_beta_range_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_beta_range (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in cb_stim_locked.
function cb_stim_locked_Callback(hObject, eventdata, handles)
% hObject    handle to cb_stim_locked (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cb_stim_locked


% --- Executes on button press in cb_broad_band.
function cb_broad_band_Callback(hObject, eventdata, handles)
% hObject    handle to cb_broad_band (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cb_broad_band



function txt_stim_freq_Callback(hObject, eventdata, handles)
% hObject    handle to txt_stim_freq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_stim_freq as text
%        str2double(get(hObject,'String')) returns contents of txt_stim_freq as a double


% --- Executes during object creation, after setting all properties.
function txt_stim_freq_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_stim_freq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in pu_metric.
function pu_metric_Callback(hObject, eventdata, handles)
% hObject    handle to pu_metric (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pu_metric contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pu_metric


% --- Executes during object creation, after setting all properties.
function pu_metric_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pu_metric (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in cb_incl_stim_freq.
function cb_incl_stim_freq_Callback(hObject, eventdata, handles)
% hObject    handle to cb_incl_stim_freq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cb_incl_stim_freq



function txt_neighb_freq_Callback(hObject, eventdata, handles)
% hObject    handle to txt_neighb_freq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_neighb_freq as text
%        str2double(get(hObject,'String')) returns contents of txt_neighb_freq as a double


% --- Executes during object creation, after setting all properties.
function txt_neighb_freq_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_neighb_freq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in cb_pred_roi.
function cb_pred_roi_Callback(hObject, eventdata, handles)
% hObject    handle to cb_pred_roi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cb_pred_roi


% --- Executes on button press in cb_roi_weight.
function cb_roi_weight_Callback(hObject, eventdata, handles)
% hObject    handle to cb_roi_weight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cb_roi_weight


% --- Executes on button press in cb_rel_sh.
function cb_rel_sh_Callback(hObject, eventdata, handles)
% hObject    handle to cb_rel_sh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cb_rel_sh


% --- Executes on button press in cb_rel_scan.
function cb_rel_scan_Callback(hObject, eventdata, handles)
% hObject    handle to cb_rel_scan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cb_rel_scan



function txt_nit_rel_Callback(hObject, eventdata, handles)
% hObject    handle to txt_nit_rel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_nit_rel as text
%        str2double(get(hObject,'String')) returns contents of txt_nit_rel as a double


% --- Executes during object creation, after setting all properties.
function txt_nit_rel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_nit_rel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txt_nit_fit_Callback(hObject, eventdata, handles)
% hObject    handle to txt_nit_fit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_nit_fit as text
%        str2double(get(hObject,'String')) returns contents of txt_nit_fit as a double


% --- Executes during object creation, after setting all properties.
function txt_nit_fit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_nit_fit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in cb_cv_roi_fit.
function cb_cv_roi_fit_Callback(hObject, eventdata, handles)
% hObject    handle to cb_cv_roi_fit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cb_cv_roi_fit

function set__widget_default_values(hObject, handles, def_params)


fnames = fieldnames(def_params);

for cur_field = fnames'
    
    switch lower(cur_field{1})
        
        case 'sigma'
            props = get(handles.pu_prf_sigma);
            set(handles.pu_prf_sigma,'Value',find_pu_val(props, def_params.(cur_field{1})));
            
        case 'sigma_fix'
            %props = get(handles.txt_sigma_fix);
            set(handles.txt_sigma_fix,'String',def_params.(cur_field{1}));
            
        case 'sigma_range'
            %props = get(handles.txt_sigma_range);
            set(handles.txt_sigma_range,'String',def_params.(cur_field{1}));
        
        case 'x0'
            props = get(handles.pu_prf_x);
            set(handles.pu_prf_x,'Value',find_pu_val(props, def_params.(cur_field{1})));
        
        case 'x0_fix'
            %props = get(handles.txt_fix_x);
            set(handles.txt_fix_x,'String',def_params.(cur_field{1}));

        case 'x0_range'
            %props = get(handles.txt_x_range);
            set(handles.txt_x_range,'String',def_params.(cur_field{1}));
            
        case 'y0'
            props = get(handles.pu_prf_y);
            set(handles.pu_prf_y,'Value',find_pu_val(props, def_params.(cur_field{1})));
            
        case 'y0_fix'
            %props = get(handles.txt_fix_y);
            set(handles.txt_fix_y,'String',def_params.(cur_field{1}));
            
        case 'y0_range'
            %props = get(handles.txt_y_range);
            set(handles.txt_y_range,'String',def_params.(cur_field{1}));

        case 'beta'
            props = get(handles.pu_prf_beta);
            set(handles.pu_prf_beta,'Value',find_pu_val(props, def_params.(cur_field{1})));
            
        case 'beta_fix'
            %props = get(handles.txt_beta_fix);
            set(handles.txt_beta_fix,'String',def_params.(cur_field{1}));
            
        case 'beta_range'
            %props = get(handles.txt_beta_range);
            set(handles.txt_beta_range,'String',def_params.(cur_field{1}));
            
        case 'roi_specific'
            %props = get(handles.cb_pred_roi);
            set(handles.cb_pred_roi,'Value',def_params.(cur_field{1}));
            
        case 'fit_roi_specific'
            %props = get(handles.cb_roi_weight);
            set(handles.cb_roi_weight,'Value',def_params.(cur_field{1}));
            
        case 'do_bb'
            %props = get(handles.cb_broad_band);
            set(handles.cb_broad_band,'Value',def_params.(cur_field{1}));
            
        case 'do_sl'
            %props = get(handles.cb_stim_locked);
            set(handles.cb_stim_locked,'Value',def_params.(cur_field{1}));
            
            
        case 'stim_freq'
            %props = get(handles.txt_stim_freq);
            set(handles.txt_stim_freq,'String',def_params.(cur_field{1}));
            
        case 'reliability_split_half'
            %props = get(handles.cb_rel_sh);
            set(handles.cb_rel_sh,'Value',def_params.(cur_field{1}));
            
        case 'reliability_scans'
            %props = get(handles.cb_rel_scan);
            set(handles.cb_rel_scan,'Value',def_params.(cur_field{1}));
            
        case 'n_iterations_rel'
            %props = get(handles.txt_nit_rel);
            set(handles.txt_nit_rel,'String',def_params.(cur_field{1}));
            
        case 'metric'
            props = get(handles.pu_metric);
            set(handles.pu_metric,'Value',find_pu_val(props, def_params.(cur_field{1})));
            
        case 'metric_include_stim_freq'
            %props = get(handles.cb_incl_stim_freq);
            set(handles.cb_incl_stim_freq,'Value',def_params.(cur_field{1}));
            
        case 'metric_nbour_range'
            %props = get(handles.txt_neighb_freq);
            set(handles.txt_neighb_freq,'String',def_params.(cur_field{1}));
            
        case 'fit_cross_validate'    
            %props = get(handles.cb_cv_roi_fit);
            set(handles.cb_cv_roi_fit,'Value',def_params.(cur_field{1}));
            
        case 'fit_n_iterations_cv'    
            %props = get(handles.txt_nit_fit);
            set(handles.txt_nit_fit,'String',def_params.(cur_field{1}));
            
        case 'comb_lr_rois'
            set(handles.cb_comb_lr_roi,'Value',def_params.(cur_field{1}));
            
            
        case 'comb_dv_rois'            
            set(handles.cb_comb_rois_dv,'Value',def_params.(cur_field{1}));
            
        case 'roi_mask'
            set(handles.cb_use_roi_mask,'Value',def_params.(cur_field{1}));
            

        case 'n_iterations_scramble'
            set(handles.txt_scramble_iterations,'String',def_params.(cur_field{1}));
            
        case 've_thr'
            set(handles.cb_ve_thr,'Value',def_params.(cur_field{1}));
            
        case 'beta_thr'
            set(handles.cb_beta_thr,'Value',def_params.(cur_field{1}));
            
        case 've_thr_vals'
            set(handles.txt_ve_thr_min,'String',def_params.(cur_field{1})(1))
            set(handles.txt_ve_thr_max,'String',def_params.(cur_field{1})(2))
        case 'beta_thr_vals'
            set(handles.txt_beta_thr_min,'String',def_params.(cur_field{1})(1))
            set(handles.txt_beta_thr_max,'String',def_params.(cur_field{1})(2))
            
            
    end
    
end


function set_val = find_pu_val(props, cur_param)

field_idx = find(strcmp(cur_param, props.String));
if isempty(field_idx)
    error('Could not find exact match for %s\n',cur_param)
else
    set_val = field_idx;
end



function data = get_val_from_widget(hObject, handles, def_params)

fnames = fieldnames(def_params);

for cur_field = fnames'
    
    switch lower(cur_field{1})
        
        case 'sigma'
            tmp = get(handles.pu_prf_sigma,'String');
            data.(cur_field{1}) = tmp{get(handles.pu_prf_sigma,'Value')};
            
        case 'sigma_fix'
            %props = get(handles.txt_sigma_fix);
            data.(cur_field{1}) = str2double(get(handles.txt_sigma_fix,'String'));
            
            
        case 'sigma_range'
            
            %props = get(handles.txt_sigma_range);
            cur_str = get(handles.txt_sigma_range,'String');
            cur_str(cur_str == '[') = '';
            cur_str(cur_str == ']') = '';
            
            if isempty(cur_str)
                data.(cur_field{1}) = [];
            else
                cel_val = strsplit(cur_str,' ');
                cel_val = cel_val(~cellfun(@isempty, cel_val));
                data.(cur_field{1}) = unique(cell2mat(cellfun(@eval,cel_val,'UniformOutput',false)));
            end
            
            
        case 'x0'
            tmp = get(handles.pu_prf_x,'String');
            data.(cur_field{1}) = tmp{get(handles.pu_prf_x,'Value')};
            
        case 'x0_fix'
            %props = get(handles.txt_fix_x);
            data.(cur_field{1}) = str2double(get(handles.txt_fix_x,'String'));
            
        case 'x0_range'
            %props = get(handles.txt_x_range);
            cur_str = get(handles.txt_x_range,'String');
            cur_str(cur_str == '[') = '';
            cur_str(cur_str == ']') = '';
            
            if isempty(cur_str)
                data.(cur_field{1}) = [];
            else
                cel_val = strsplit(cur_str,' ');
                cel_val = cel_val(~cellfun(@isempty, cel_val));
                data.(cur_field{1}) = unique(cell2mat(cellfun(@eval,cel_val,'UniformOutput',false)));
            end
            
            
            
        case 'y0'
            tmp = get(handles.pu_prf_y,'String');
            data.(cur_field{1}) = tmp{get(handles.pu_prf_y,'Value')};
            
        case 'y0_fix'
            %props = get(handles.txt_fix_y);
            data.(cur_field{1}) = str2double(get(handles.txt_fix_y,'String'));
            
        case 'y0_range'
            %props = get(handles.txt_y_range);
            cur_str = get(handles.txt_y_range,'String');
            cur_str(cur_str == '[') = '';
            cur_str(cur_str == ']') = '';
            
            if isempty(cur_str)
                data.(cur_field{1}) = [];
            else
                cel_val = strsplit(cur_str,' ');
                cel_val = cel_val(~cellfun(@isempty, cel_val));
                data.(cur_field{1}) = unique(cell2mat(cellfun(@eval,cel_val,'UniformOutput',false)));
            end
            
            
        case 'beta'
            tmp = get(handles.pu_prf_beta,'String');
            data.(cur_field{1}) = tmp{get(handles.pu_prf_beta,'Value')};            
            
        case 'beta_fix'
            %props = get(handles.txt_beta_fix);
            data.(cur_field{1}) = str2double(get(handles.txt_beta_fix,'String'));
            
        case 'beta_range'
            %props = get(handles.txt_beta_range);
            cur_str = get(handles.txt_beta_range,'String');
            cur_str(cur_str == '[') = '';
            cur_str(cur_str == ']') = '';
            
            
            if isempty(cur_str)
                data.(cur_field{1}) = [];
            else
                cel_val = strsplit(cur_str,' ');
                cel_val = cel_val(~cellfun(@isempty, cel_val));
                data.(cur_field{1}) = unique(cell2mat(cellfun(@eval,cel_val,'UniformOutput',false)));
            end
            
        case 'roi_specific'
            %props = get(handles.cb_pred_roi);            
            data.(cur_field{1}) = get(handles.cb_pred_roi,'Value');            
            
        case 'fit_roi_specific'
            %props = get(handles.cb_roi_weight);
            data.(cur_field{1}) = get(handles.cb_roi_weight,'Value');            
            
            
        case 'do_bb'
            %props = get(handles.cb_broad_band);
            data.(cur_field{1}) = get(handles.cb_broad_band,'Value');            
            
        case 'do_sl'
            %props = get(handles.cb_stim_locked);
            data.(cur_field{1}) = get(handles.cb_stim_locked,'Value');            
            
        case 'stim_freq'
            %props = get(handles.txt_stim_freq);
            data.(cur_field{1}) = str2double(get(handles.txt_stim_freq,'String'));
            
        case 'reliability_split_half'
            %props = get(handles.cb_rel_sh);
            data.(cur_field{1}) = get(handles.cb_rel_sh,'Value');            
            
            
        case 'reliability_scans'
            %props = get(handles.cb_rel_scan);
            data.(cur_field{1}) = get(handles.cb_rel_scan,'Value');
            
        case 'n_iterations_rel'
            %props = get(handles.txt_nit_rel);
            data.(cur_field{1}) = str2double(get(handles.txt_nit_rel,'String'));
            
        case 'metric'
            tmp = get(handles.pu_metric,'String');
            data.(cur_field{1}) = tmp{get(handles.pu_metric,'Value')};
            
            
        case 'metric_include_stim_freq'
            %props = get(handles.cb_incl_stim_freq);
            data.(cur_field{1}) = get(handles.cb_incl_stim_freq,'Value');
            
        case 'metric_nbour_range'
            %props = get(handles.txt_neighb_freq);
            cur_str = get(handles.txt_neighb_freq,'String');
            cur_str(cur_str == '[') = '';
            cur_str(cur_str == ']') = '';

            
            if isempty(cur_str)
                data.(cur_field{1}) = [];
            else
                cel_val = strsplit(cur_str,' ');
                cel_val = cel_val(~cellfun(@isempty, cel_val));
                data.(cur_field{1}) = unique(cell2mat(cellfun(@eval,cel_val,'UniformOutput',false)));
            end
            
            
        case 'fit_cross_validate'    
            %props = get(handles.cb_cv_roi_fit);
            data.(cur_field{1}) = get(handles.cb_cv_roi_fit,'Value');
            
            
        case 'fit_n_iterations_cv'
            %props = get(handles.txt_nit_fit);
            data.(cur_field{1}) = str2double(get(handles.txt_nit_fit,'String'));
            
            
        case 'comb_lr_rois'    
            data.(cur_field{1}) = get(handles.cb_comb_lr_roi,'Value');
            
        case 'comb_dv_rois'
            data.(cur_field{1}) = get(handles.cb_comb_rois_dv,'Value');
            
        case 'roi_mask'
            data.(cur_field{1}) = get(handles.cb_use_roi_mask,'Value');
            
        case 'n_iterations_scramble'
            data.(cur_field{1}) = str2double(get(handles.txt_scramble_iterations,'String'));
            
            
        case 've_thr'
            data.(cur_field{1}) = get(handles.cb_ve_thr,'Value');
            
        case 'beta_thr'
            data.(cur_field{1}) = get(handles.cb_beta_thr,'Value');
            
        case 've_thr_vals'
            data.(cur_field{1}) = [str2double(get(handles.txt_ve_thr_min,'String')) ...
                str2double(get(handles.txt_ve_thr_max,'String'))];
            
        case 'beta_thr_vals'
            data.(cur_field{1}) = [str2double(get(handles.txt_beta_thr_min,'String')) ...
                str2double(get(handles.txt_beta_thr_max,'String'))];
            
            
            
            
            
    end
    
end

function have_all = check_parameter_out_fields(handles)
have_all = false;
missing = cell(1,20);
params_out = handles.ch_model.params_out;

cnt = 0;

if ~isempty(params_out.x0) && ischar(params_out.x0)
else
    cnt = cnt+1;
    missing{cnt} = 'pRF X0';
end

if ~isempty(params_out.y0) && ischar(params_out.y0)
else
    cnt = cnt+1;
    missing{cnt} = 'pRF y0';
end

if ~isempty(params_out.sigma) && ischar(params_out.sigma)
else
    cnt = cnt+1;
    missing{cnt} = 'pRF sigma';
end

if ~isempty(params_out.beta) && ischar(params_out.beta)
else
    cnt = cnt+1;
    missing{cnt} = 'pRF beta';
end

if strcmpi(params_out.x0, 'fixed x (absolute)') || ...
        strcmpi(params_out.x0, 'fixed x (proportional)') ...
        strcmpi(params_out.x0, 'fixed x (proportional, smoothed)')
    
    if ~isempty(params_out.x0_fix) && ~any(isnan(params_out.x0_fix))
        
    else
        cnt = cnt+1;
        missing{cnt} = 'Fixed X0';
    end
    
    
end

if strcmpi(params_out.y0, 'fixed y (absolute)') || ...
        strcmpi(params_out.y0, 'fixed y (proportion)') ...
        strcmpi(params_out.y0, 'fixed y (proportion, smoothed)')
    
    if ~isempty(params_out.y0_fix) && ~any(isnan(params_out.y0_fix))
        
    else
        cnt = cnt+1;
        missing{cnt} = 'Fixed Y0';
    end
    
    
end

if strcmpi(params_out.sigma, 'fixed sigma (absolute)') || ...
        strcmpi(params_out.sigma, 'fixed sigma (proportion)') ...
        strcmpi(params_out.sigma, 'fixed sigma (proportion, smoothed)')
    
    if ~isempty(params_out.sigma_fix) && ~any(isnan(params_out.sigma_fix))
        
    else
        cnt = cnt+1;
        missing{cnt} = 'Fixed Sigma';
    end
    
    
end

if strcmpi(params_out.beta, 'fixed mresp_smoothed (absolute)') || ...
    strcmpi(params_out.beta, 'fixed mresp_smoothed (proportion)') || ... 
    strcmpi(params_out.beta, 'fixed mresp (absolute)') || ...
    strcmpi(params_out.beta, 'fixed mresp (proportion)') || ...
    strcmpi(params_out.beta, 'fixed recomp_beta (absolute)') || ...
    strcmpi(params_out.beta, 'fixed recomp_beta (proportion)') || ...
    strcmpi(params_out.beta, 'fixed beta (absolute)') || ...
    strcmpi(params_out.beta, 'fixed beta (proportion)')

        
    if ~isempty(params_out.beta_fix) && ~any(isnan(params_out.beta_fix))
        
    else
        cnt = cnt+1;
        missing{cnt} = 'Fixed Beta';
    end

end



if strcmpi(params_out.x0, 'x range(absolute)') || ...
        strcmpi(params_out.x0, 'x range(proportional)') ...
        strcmpi(params_out.x0, 'x range(proportional, smoothed)')
    
    
    if ~isempty(params_out.x0_range) && ~any(isnan(params_out.x0_range))
        
    else
        cnt = cnt+1;
        missing{cnt} = 'x range';
    end
    
end

if strcmpi(params_out.y0, 'y range (absolute)') || ...
        strcmpi(params_out.y0, 'y range (proportion)') ...
        strcmpi(params_out.y0, 'y range (proportion, smoothed)')
    
    if ~isempty(params_out.y0_range) && ~any(isnan(params_out.y0_range))
        
    else
        cnt = cnt+1;
        missing{cnt} = 'y range';
    end
    
end

if strcmpi(params_out.sigma, 'sigma range (absolute)') || ...
        strcmpi(params_out.sigma, 'sigma range (proportion)') ...
        strcmpi(params_out.sigma, 'sigma range (proportion, smoothed)')
    
    if ~isempty(params_out.sigma_range) && ~any(isnan(params_out.sigma_range))
        
    else
        cnt = cnt+1;
        missing{cnt} = 'Sigma range';
    end
    
end

if strcmpi(params_out.beta, 'mresp_smoothed range (absolute)') || ...
        strcmpi(params_out.beta, 'mresp_smoothed range (proportion)') || ...
        strcmpi(params_out.beta, 'mresp range (absolute)') || ...
        strcmpi(params_out.beta, 'mresp range (proportion)') || ...
        strcmpi(params_out.beta, 'recomp_beta range (absolute)') || ...
        strcmpi(params_out.beta, 'recomp_beta range (proportion)') || ...
        strcmpi(params_out.beta, 'beta range (absolute)') || ...
        strcmpi(params_out.beta, 'beta range (proportion)')
    
    if ~isempty(params_out.beta_range) && ~any(isnan(params_out.beta_range))
        
    else
        cnt = cnt+1;
        missing{cnt} = 'Beta range';
    end
    
    
    
end

if any(strfind(params_out.sigma, 'scrambled')) || ...
        any(strfind(params_out.x0, 'scrambled')) || ...
        any(strfind(params_out.y0, 'scrambled')) || ...
        any(strfind(params_out.beta, 'scrambled'))
    
    if ~isempty(params_out.n_iterations_scramble) ..., 
            && ~any(isnan(params_out.metric_nbour_range))
        
    else
        cnt = cnt+1;
        missing{cnt} = 'Scramble iterations';
    end
    
    
end



if params_out.do_sl || params_out.do_bb
else
    cnt = cnt+1;
    missing{cnt} = 'Stimulus locked or broad band, or both';
end

if strcmp(params_out.metric,'stimulus/neighbour frequencies')
    if ~isempty(params_out.metric_nbour_range) && ~any(isnan(params_out.metric_nbour_range))
                
    else
        cnt = cnt+1;
            missing{cnt} = 'Neighbouring frequencies for metric';
    end
    
end

if ~isempty(params_out.stim_freq) && ~isnan(params_out.stim_freq) ...
        && ~isempty(params_out.stim_freq)
else
    cnt = cnt+1;
    missing{cnt} = 'Stimulus frequency';
end


if params_out.fit_cross_validate
    if ~isempty(params_out.fit_n_iterations_cv) && ~isnan(params_out.fit_n_iterations_cv)
        
    else
        cnt = cnt+1;
        missing{cnt} = 'N iterations cross validation';
    end
end


if params_out.reliability_scans || params_out.reliability_split_half
    if ~isempty(params_out.n_iterations_rel) && ~isnan(params_out.n_iterations_rel)
        
    else
        cnt = cnt+1;
        missing{cnt} = 'N iterations reliability';
    end
end


switch lower(handles.ch_model.type)


    case 'run original model'

    case 'equal weight'
        if ~isempty(params_out.beta_fix) && ~isnan(params_out.beta_fix) ...
                && ~isempty(params_out.beta_fix)
        else
            cnt = cnt+1;
            missing{cnt} = 'pRF beta fix';
        end
    case 'scramble prf parameters'
           
    case 'fix prf size'
        
        if ~isempty(params_out.sigma_fix) && ~isnan(params_out.sigma_fix) ...
                && ~isempty(params_out.sigma_fix)
        else
            cnt = cnt+1;
            missing{cnt} = 'pRF sigma fix';
        end
        
        
    case 'prf size range'
        
        if ~isempty(params_out.sigma_range) && ~isnan(params_out.sigma_range) ...
                && ~isempty(params_out.sigma_range)
        else
            cnt = cnt+1;
            missing{cnt} = 'pRF sigma range';
        end
        
        
           
    case 'reliability check'
                
        if params_out.reliability_scans || params_out.reliability_split_half
        else
            cnt = cnt+1;
            missing{cnt} = 'Reliability split or scan, or both';
        end
         
        
        if ~isempty(params_out.n_iterations_rel) && ~isnan(params_out.n_iterations_rel)
        else
            cnt = cnt+1;
            missing{cnt} = 'N iterations reliability';
        end
        
        
    case 'fit separate roi predictions'
        if params_out.roi_specific
        else
            cnt = cnt+1;
            missing{cnt} = 'Predictions per ROI';
        end
        
        if params_out.fit_roi_specific	
        else
            cnt = cnt+1;
            missing{cnt} = 'Fit predictions per ROI';
        end
        
end

if ~any(~cellfun(@isempty, missing))
   have_all = true;
else
   missing = missing(1:cnt); 
   
end

if have_all
    
else
    fprintf('Please provide valid values for:\n')
    for n = 1:length(missing)
        fprintf('%s\n',missing{n})
    end
        
    
end


% --- Executes on button press in pb_done_me.
function pb_done_me_Callback(hObject, eventdata, handles)
% hObject    handle to pb_done_me (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.ch_model.params_out = get_val_from_widget(hObject, handles, handles.ch_model.def_params);

go_on = check_parameter_out_fields(handles);
guidata(hObject, handles);
if go_on
    figure1_CloseRequestFcn(handles.figure1, eventdata, handles);
end

% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles) %#ok<*INUSD>
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isequal(get(hObject, 'waitstatus'), 'waiting')
    uiresume(hObject);
    
else
    % Hint: delete(hObject) closes the figure
    delete(hObject);
end


% --- Executes when selected object is changed in bg_radio_beta.
function bg_radio_beta_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in bg_radio_beta 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

get(eventdata.NewValue,'String')


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over rb_equal_beta_weight.
function rb_equal_beta_weight_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to rb_equal_beta_weight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in cb_use_roi_mask.
function cb_use_roi_mask_Callback(hObject, eventdata, handles)
% hObject    handle to cb_use_roi_mask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cb_use_roi_mask


% --- Executes on button press in cb_comb_lr_roi.
function cb_comb_lr_roi_Callback(hObject, eventdata, handles)
% hObject    handle to cb_comb_lr_roi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cb_comb_lr_roi


% --- Executes on button press in cb_comb_rois_dv.
function cb_comb_rois_dv_Callback(hObject, eventdata, handles)
% hObject    handle to cb_comb_rois_dv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cb_comb_rois_dv



function edit13_Callback(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit13 as text
%        str2double(get(hObject,'String')) returns contents of edit13 as a double


% --- Executes during object creation, after setting all properties.
function edit13_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txt_scramble_iterations_Callback(hObject, eventdata, handles)
% hObject    handle to txt_scramble_iterations (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_scramble_iterations as text
%        str2double(get(hObject,'String')) returns contents of txt_scramble_iterations as a double


% --- Executes during object creation, after setting all properties.
function txt_scramble_iterations_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_scramble_iterations (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in cb_ve_thr.
function cb_ve_thr_Callback(hObject, eventdata, handles)
% hObject    handle to cb_ve_thr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cb_ve_thr


% --- Executes on button press in cb_beta_thr.
function cb_beta_thr_Callback(hObject, eventdata, handles)
% hObject    handle to cb_beta_thr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cb_beta_thr



function txt_ve_thr_min_Callback(hObject, eventdata, handles)
% hObject    handle to txt_ve_thr_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_ve_thr_min as text
%        str2double(get(hObject,'String')) returns contents of txt_ve_thr_min as a double


% --- Executes during object creation, after setting all properties.
function txt_ve_thr_min_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_ve_thr_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txt_ve_thr_max_Callback(hObject, eventdata, handles)
% hObject    handle to txt_ve_thr_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_ve_thr_max as text
%        str2double(get(hObject,'String')) returns contents of txt_ve_thr_max as a double


% --- Executes during object creation, after setting all properties.
function txt_ve_thr_max_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_ve_thr_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txt_beta_thr_max_Callback(hObject, eventdata, handles)
% hObject    handle to txt_beta_thr_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_beta_thr_max as text
%        str2double(get(hObject,'String')) returns contents of txt_beta_thr_max as a double


% --- Executes during object creation, after setting all properties.
function txt_beta_thr_max_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_beta_thr_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txt_beta_thr_min_Callback(hObject, eventdata, handles)
% hObject    handle to txt_beta_thr_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_beta_thr_min as text
%        str2double(get(hObject,'String')) returns contents of txt_beta_thr_min as a double


% --- Executes during object creation, after setting all properties.
function txt_beta_thr_min_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_beta_thr_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
