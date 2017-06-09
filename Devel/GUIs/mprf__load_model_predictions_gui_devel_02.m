function varargout = mprf__load_model_predictions_gui_devel_02(varargin)
% MPRF__LOAD_MODEL_PREDICTIONS_GUI_DEVEL_02 MATLAB code for mprf__load_model_predictions_gui_devel_02.fig
%      MPRF__LOAD_MODEL_PREDICTIONS_GUI_DEVEL_02, by itself, creates a new MPRF__LOAD_MODEL_PREDICTIONS_GUI_DEVEL_02 or raises the existing
%      singleton*.
%
%      H = MPRF__LOAD_MODEL_PREDICTIONS_GUI_DEVEL_02 returns the handle to a new MPRF__LOAD_MODEL_PREDICTIONS_GUI_DEVEL_02 or the handle to
%      the existing singleton*.
%
%      MPRF__LOAD_MODEL_PREDICTIONS_GUI_DEVEL_02('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MPRF__LOAD_MODEL_PREDICTIONS_GUI_DEVEL_02.M with the given input arguments.
%
%      MPRF__LOAD_MODEL_PREDICTIONS_GUI_DEVEL_02('Property','Value',...) creates a new MPRF__LOAD_MODEL_PREDICTIONS_GUI_DEVEL_02 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before mprf__load_model_predictions_gui_devel_02_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to mprf__load_model_predictions_gui_devel_02_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help mprf__load_model_predictions_gui_devel_02

% Last Modified by GUIDE v2.5 09-Jun-2017 14:39:28
% ADD:
% IF ROI SPECIFIC OR NOT
% AMOUNT OF ROIS
% STIMULUS WINDOW
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @mprf__load_model_predictions_gui_devel_02_OpeningFcn, ...
                   'gui_OutputFcn',  @mprf__load_model_predictions_gui_devel_02_OutputFcn, ...
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


% --- Executes just before mprf__load_model_predictions_gui_devel_02 is made visible.
function mprf__load_model_predictions_gui_devel_02_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to mprf__load_model_predictions_gui_devel_02 (see VARARGIN)

% Choose default command line output for mprf__load_model_predictions_gui_devel_02
handles.output = hObject;

handles.dir.main = mprf__get_directory('main_dir');
handles.dir.model_pred = mprf__get_directory('model_predictions');
tmp = dir(fullfile(handles.dir.main, handles.dir.model_pred, '*.mat'));
handles.models.names = {tmp.name};

set(handles.lst_models_available, 'String',handles.models.names)

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes mprf__load_model_predictions_gui_devel_02 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = mprf__load_model_predictions_gui_devel_02_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in lst_models_available.
function lst_models_available_Callback(hObject, eventdata, handles)
% hObject    handle to lst_models_available (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns lst_models_available contents as cell array
%        contents{get(hObject,'Value')} returns selected item from lst_models_available

handles.selected.model = fullfile(handles.dir.main, ...
    handles.dir.model_pred, ... 
    handles.models.names{get(hObject,'Value')});


handles.tmp.model = load(handles.selected.model,'model','roi','prf');

beta_opts = handles.tmp.model.prf.beta;
x0_opts = handles.tmp.model.prf.x0;
y0_opts = handles.tmp.model.prf.y0;
sigma_opts = handles.tmp.model.prf.sigma;


set(handles.txt_x0_par,'String',handles.tmp.model.model.params.x0);
set(handles.txt_y0_par,'String',handles.tmp.model.model.params.y0);
set(handles.txt_sigma_par,'String',handles.tmp.model.model.params.sigma);
set(handles.txt_beta_par,'String',handles.tmp.model.model.params.beta);

set(handles.txt_beta_range,'String',get_prf_par_val(handles, beta_opts, 'beta'));
set(handles.txt_y0_range,'String',get_prf_par_val(handles, y0_opts, 'y0'));
set(handles.txt_x0_range,'String',get_prf_par_val(handles, x0_opts, 'x0'));
set(handles.txt_sigma_range,'String',get_prf_par_val(handles, sigma_opts, 'sigma'));



guidata(hObject, handles)

% --- Executes during object creation, after setting all properties.
function lst_models_available_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lst_models_available (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function val = get_prf_par_val(handles, opts,par_name)

if opts.type.fixed
    val = num2str(handles.tmp.model.model.params.([par_name '_fix']));
    
elseif opts.type.range
    tmp = handles.tmp.model.model.params.([par_name '_range']);
    val = [num2str(min(tmp)) ':' num2str(diff(tmp([1 2]))) ':' num2str(max(tmp))];
    
else
    val = 'NA';
    
end
