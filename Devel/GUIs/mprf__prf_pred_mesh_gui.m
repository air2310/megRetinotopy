function varargout = mprf__prf_pred_mesh_gui(varargin)
% MPRF__PRF_PRED_MESH_GUI MATLAB code for mprf__prf_pred_mesh_gui.fig
%      MPRF__PRF_PRED_MESH_GUI, by itself, creates a new MPRF__PRF_PRED_MESH_GUI or raises the existing
%      singleton*.
%
%      H = MPRF__PRF_PRED_MESH_GUI returns the handle to a new MPRF__PRF_PRED_MESH_GUI or the handle to
%      the existing singleton*.
%
%      MPRF__PRF_PRED_MESH_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MPRF__PRF_PRED_MESH_GUI.M with the given input arguments.
%
%      MPRF__PRF_PRED_MESH_GUI('Property','Value',...) creates a new MPRF__PRF_PRED_MESH_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before mprf__prf_pred_mesh_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to mprf__prf_pred_mesh_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help mprf__prf_pred_mesh_gui

% Last Modified by GUIDE v2.5 05-Jul-2017 14:58:15

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @mprf__prf_pred_mesh_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @mprf__prf_pred_mesh_gui_OutputFcn, ...
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


% --- Executes just before mprf__prf_pred_mesh_gui is made visible.
function mprf__prf_pred_mesh_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to mprf__prf_pred_mesh_gui (see VARARGIN)

handles.data = varargin{1};
handles.mesh.surf_file = varargin{2};
handles.data.stimulus.full_im = handles.data.stimulus.full_im./max(handles.data.stimulus.full_im(:));
fnames = fieldnames(handles);
cb_nums = nan(size(fnames));
ctr = 0;

for n = 1:length(fnames)
   cur_fname = fnames{n};
   if any(strfind(cur_fname, 'roi_checkbox')) 
       ctr = ctr+1;
       cb_nums(ctr) = str2double(cur_fname(length('roi_checkbox')+1:end));
   end    
end
cb_nums = cb_nums(1:ctr);
cb_nums = sort(cb_nums);

if isfield(handles.data.roi,'idx')
   roi_names = fieldnames(handles.data.roi.idx);
    
   for nn = 1:length(roi_names)
      eval(['set(handles.roi_checkbox' num2str(cb_nums(nn)) ',''String'',roi_names{nn});']) 
       
   end
    
else
    nn=1;
    
    eval(['set(handles.roi_checkbox' num2str(cb_nums(nn)) ',''String'',''All'');'])
    
end

for nnn = nn+1:length(cb_nums)
    eval(['set(handles.roi_checkbox' num2str(cb_nums(nnn)) ',''String'',''Empty'',''Enable'',''off'');'])
    
    
end

eval(['set(handles.roi_checkbox' num2str(cb_nums(1)) ',''Value'',1);'])
stim_range = [1 size(handles.data.stimulus.full_im,3)];

set(handles.stim_slider,...
    'Min',stim_range(1),...
    'Max',stim_range(2),...
    'SliderStep',  [(1/(stim_range(2) - stim_range(1))) (10/(stim_range(2) - stim_range(1)))],...
    'Value',1);

set(handles.txt_cur_stim,'String','1');
set(handles.txt_max_stim,'String',num2str(stim_range(2)));
axes(handles.axes1)

colormap('gray');

handles = updateStim(hObject, handles);
handles = initMesh(handles);
handles.cb_nums = cb_nums(1:nn);
guidata(hObject, handles);

updateMesh(hObject, handles);

% Choose default command line output for mprf__prf_pred_mesh_gui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes mprf__prf_pred_mesh_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


function updateMesh(hObject, handles)

% Get stim position
cur_im = ceil(get(handles.stim_slider,'Value'));
% Get ROIs
cur_cbs = handles.cb_nums;

[n_im, n_vert, n_roi] = size(handles.data.pred_resp);

all_rois = zeros(size(cur_cbs));
ctr = 0;

pred = zeros(1,n_vert);
mask = nan(size(pred));

for nn = 1:length(cur_cbs)
    eval(['cur_val = get(handles.roi_checkbox' num2str(cur_cbs(nn)) ',''Value'');'])
    
    if cur_val
        ctr = ctr+1;
        eval(['cur_roi_str = get(handles.roi_checkbox' num2str(cur_cbs(nn)) ',''String'');'])
        all_rois(ctr) = nn;
        
        pred = sum([pred; squeeze(handles.data.pred_resp(cur_im,:,nn))]);
        mask(handles.data.roi.idx_out == handles.data.roi.idx.(cur_roi_str).idx) = true;
        
    end
    
end
drange = [];
handles.mesh.bs_msh = mprfSessionColorMesh(handles.mesh.bs_msh, pred,jet(256), drange, ~isnan(mask));




guidata(hObject, handles);




function handles = initMesh(handles)

handles.mesh.bs_msh = mprfMeshFromBrainstorm(handles.mesh.surf_file);
handles.mesh.bs_msh = meshVisualize(handles.mesh.bs_msh);



function handles = updateStim(hObject, handles)

cur_val = ceil(get(handles.stim_slider,'Value'));

axes(handles.axes1);
imagesc(handles.data.stimulus.full_im(:,:,cur_val), [0 1]);
set(handles.axes1,'YTick',[]);
set(handles.axes1,'XTick',[]);
set(handles.txt_cur_stim, 'String',num2str(cur_val));







% --- Outputs from this function are returned to the command line.
function varargout = mprf__prf_pred_mesh_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

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

updateStim(hObject, handles);
updateMesh(hObject, handles);

% --- Executes during object creation, after setting all properties.
function stim_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to stim_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in roi_checkbox22.
function roi_checkbox22_Callback(hObject, eventdata, handles)
% hObject    handle to roi_checkbox22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of roi_checkbox22
updateMesh(hObject, handles)

% --- Executes on button press in roi_checkbox23.
function roi_checkbox23_Callback(hObject, eventdata, handles)
% hObject    handle to roi_checkbox23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of roi_checkbox23
updateMesh(hObject, handles)

% --- Executes on button press in roi_checkbox24.
function roi_checkbox24_Callback(hObject, eventdata, handles)
% hObject    handle to roi_checkbox24 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of roi_checkbox24
updateMesh(hObject, handles)

% --- Executes on button press in roi_checkbox25.
function roi_checkbox25_Callback(hObject, eventdata, handles)
% hObject    handle to roi_checkbox25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of roi_checkbox25
updateMesh(hObject, handles)

% --- Executes on button press in roi_checkbox26.
function roi_checkbox26_Callback(hObject, eventdata, handles)
% hObject    handle to roi_checkbox26 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of roi_checkbox26
updateMesh(hObject, handles)

% --- Executes on button press in roi_checkbox27.
function roi_checkbox27_Callback(hObject, eventdata, handles)
% hObject    handle to roi_checkbox27 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of roi_checkbox27
updateMesh(hObject, handles)

% --- Executes on button press in roi_checkbox28.
function roi_checkbox28_Callback(hObject, eventdata, handles)
% hObject    handle to roi_checkbox28 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of roi_checkbox28
updateMesh(hObject, handles)

% --- Executes on button press in roi_checkbox29.
function roi_checkbox29_Callback(hObject, eventdata, handles)
% hObject    handle to roi_checkbox29 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of roi_checkbox29
updateMesh(hObject, handles)

% --- Executes on button press in roi_checkbox30.
function roi_checkbox30_Callback(hObject, eventdata, handles)
% hObject    handle to roi_checkbox30 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of roi_checkbox30
updateMesh(hObject, handles)

% --- Executes on button press in roi_checkbox31.
function roi_checkbox31_Callback(hObject, eventdata, handles)
% hObject    handle to roi_checkbox31 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of roi_checkbox31
updateMesh(hObject, handles)

% --- Executes on button press in roi_checkbox32.
function roi_checkbox32_Callback(hObject, eventdata, handles)
% hObject    handle to roi_checkbox32 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of roi_checkbox32
updateMesh(hObject, handles)

% --- Executes on button press in roi_checkbox33.
function roi_checkbox33_Callback(hObject, eventdata, handles)
% hObject    handle to roi_checkbox33 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of roi_checkbox33
updateMesh(hObject, handles)

% --- Executes on button press in roi_checkbox34.
function roi_checkbox34_Callback(hObject, eventdata, handles)
% hObject    handle to roi_checkbox34 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of roi_checkbox34
updateMesh(hObject, handles)

% --- Executes on button press in roi_checkbox35.
function roi_checkbox35_Callback(hObject, eventdata, handles)
% hObject    handle to roi_checkbox35 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of roi_checkbox35
updateMesh(hObject, handles)

% --- Executes on button press in roi_checkbox36.
function roi_checkbox36_Callback(hObject, eventdata, handles)
% hObject    handle to roi_checkbox36 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of roi_checkbox36
updateMesh(hObject, handles)

% --- Executes on button press in roi_checkbox37.
function roi_checkbox37_Callback(hObject, eventdata, handles)
% hObject    handle to roi_checkbox37 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of roi_checkbox37
updateMesh(hObject, handles)

% --- Executes on button press in roi_checkbox38.
function roi_checkbox38_Callback(hObject, eventdata, handles)
% hObject    handle to roi_checkbox38 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of roi_checkbox38
updateMesh(hObject, handles)

% --- Executes on button press in roi_checkbox39.
function roi_checkbox39_Callback(hObject, eventdata, handles)
% hObject    handle to roi_checkbox39 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of roi_checkbox39
updateMesh(hObject, handles)

% --- Executes on button press in roi_checkbox40.
function roi_checkbox40_Callback(hObject, eventdata, handles)
% hObject    handle to roi_checkbox40 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of roi_checkbox40
updateMesh(hObject, handles)

% --- Executes on button press in roi_checkbox41.
function roi_checkbox41_Callback(hObject, eventdata, handles)
% hObject    handle to roi_checkbox41 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of roi_checkbox41
updateMesh(hObject, handles)

% --- Executes on button press in roi_checkbox42.
function roi_checkbox42_Callback(hObject, eventdata, handles)
% hObject    handle to roi_checkbox42 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of roi_checkbox42
updateMesh(hObject, handles)

% --- Executes on button press in roi_checkbox43.
function roi_checkbox43_Callback(hObject, eventdata, handles)
% hObject    handle to roi_checkbox43 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of roi_checkbox43
updateMesh(hObject, handles)

% --- Executes on button press in roi_checkbox44.
function roi_checkbox44_Callback(hObject, eventdata, handles)
% hObject    handle to roi_checkbox44 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of roi_checkbox44
updateMesh(hObject, handles)

% --- Executes on button press in roi_checkbox45.
function roi_checkbox45_Callback(hObject, eventdata, handles)
% hObject    handle to roi_checkbox45 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of roi_checkbox45
updateMesh(hObject, handles)

% --- Executes on button press in roi_checkbox46.
function roi_checkbox46_Callback(hObject, eventdata, handles)
% hObject    handle to roi_checkbox46 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of roi_checkbox46
updateMesh(hObject, handles)

% --- Executes on button press in roi_checkbox47.
function roi_checkbox47_Callback(hObject, eventdata, handles)
% hObject    handle to roi_checkbox47 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of roi_checkbox47
updateMesh(hObject, handles)

% --- Executes on button press in roi_checkbox48.
function roi_checkbox48_Callback(hObject, eventdata, handles)
% hObject    handle to roi_checkbox48 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of roi_checkbox48
updateMesh(hObject, handles)

% --- Executes on button press in roi_checkbox49.
function roi_checkbox49_Callback(hObject, eventdata, handles)
% hObject    handle to roi_checkbox49 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of roi_checkbox49
updateMesh(hObject, handles)

% --- Executes on button press in roi_checkbox50.
function roi_checkbox50_Callback(hObject, eventdata, handles)
% hObject    handle to roi_checkbox50 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of roi_checkbox50
updateMesh(hObject, handles)

% --- Executes on button press in roi_checkbox51.
function roi_checkbox51_Callback(hObject, eventdata, handles)
% hObject    handle to roi_checkbox51 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of roi_checkbox51
updateMesh(hObject, handles)

% --- Executes on button press in roi_checkbox52.
function roi_checkbox52_Callback(hObject, eventdata, handles)
% hObject    handle to roi_checkbox52 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of roi_checkbox52
updateMesh(hObject, handles)

% --- Executes on button press in roi_checkbox53.
function roi_checkbox53_Callback(hObject, eventdata, handles)
% hObject    handle to roi_checkbox53 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of roi_checkbox53
updateMesh(hObject, handles)

% --- Executes on button press in roi_checkbox54.
function roi_checkbox54_Callback(hObject, eventdata, handles)
% hObject    handle to roi_checkbox54 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of roi_checkbox54
updateMesh(hObject, handles)

% --- Executes on button press in roi_checkbox55.
function roi_checkbox55_Callback(hObject, eventdata, handles)
% hObject    handle to roi_checkbox55 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of roi_checkbox55
updateMesh(hObject, handles)

% --- Executes on button press in roi_checkbox56.
function roi_checkbox56_Callback(hObject, eventdata, handles)
% hObject    handle to roi_checkbox56 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of roi_checkbox56
updateMesh(hObject, handles)

% --- Executes on button press in roi_checkbox57.
function roi_checkbox57_Callback(hObject, eventdata, handles)
% hObject    handle to roi_checkbox57 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of roi_checkbox57
updateMesh(hObject, handles)

% --- Executes on button press in roi_checkbox58.
function roi_checkbox58_Callback(hObject, eventdata, handles)
% hObject    handle to roi_checkbox58 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of roi_checkbox58
updateMesh(hObject, handles)

% --- Executes on button press in roi_checkbox59.
function roi_checkbox59_Callback(hObject, eventdata, handles)
% hObject    handle to roi_checkbox59 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of roi_checkbox59
updateMesh(hObject, handles)

% --- Executes on button press in roi_checkbox60.
function roi_checkbox60_Callback(hObject, eventdata, handles)
% hObject    handle to roi_checkbox60 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of roi_checkbox60
updateMesh(hObject, handles)

% --- Executes on button press in roi_checkbox61.
function roi_checkbox61_Callback(hObject, eventdata, handles)
% hObject    handle to roi_checkbox61 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of roi_checkbox61
updateMesh(hObject, handles)