function varargout = mprf__load_model_predictions_gui(varargin)
% MPRF__LOAD_MODEL_PREDICTIONS_GUI MATLAB code for mprf__load_model_predictions_gui.fig
%      MPRF__LOAD_MODEL_PREDICTIONS_GUI, by itself, creates a new MPRF__LOAD_MODEL_PREDICTIONS_GUI or raises the existing
%      singleton*.
%
%      H = MPRF__LOAD_MODEL_PREDICTIONS_GUI returns the handle to a new MPRF__LOAD_MODEL_PREDICTIONS_GUI or the handle to
%      the existing singleton*.
%
%      MPRF__LOAD_MODEL_PREDICTIONS_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MPRF__LOAD_MODEL_PREDICTIONS_GUI.M with the given input arguments.
%
%      MPRF__LOAD_MODEL_PREDICTIONS_GUI('Property','Value',...) creates a new MPRF__LOAD_MODEL_PREDICTIONS_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before mprf__load_model_predictions_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to mprf__load_model_predictions_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help mprf__load_model_predictions_gui

% Last Modified by GUIDE v2.5 08-Jun-2017 14:30:12

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @mprf__load_model_predictions_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @mprf__load_model_predictions_gui_OutputFcn, ...
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


% --- Executes just before mprf__load_model_predictions_gui is made visible.
function mprf__load_model_predictions_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to mprf__load_model_predictions_gui (see VARARGIN)

% Choose default command line output for mprf__load_model_predictions_gui
handles.output = hObject;

handles.dir.main = mprf__get_directory('main_dir');
handles.dir.model_pred = mprf__get_directory('model_predictions');
tmp = dir(fullfile(handles.dir.main, handles.dir.model_pred, '*.mat'));
handles.models.names = {tmp.name};

set(handles.lst_models_available, 'String',handles.models.names)

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes mprf__load_model_predictions_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = mprf__load_model_predictions_gui_OutputFcn(hObject, eventdata, handles) 
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
    cellstr(get(hObject,'String')));


handles.tmp.model = load(handles.selected.model{1},'model','roi');
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


