function varargout = mbc_vmfmm(varargin)

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @mbc_vmfmm_OpeningFcn, ...
                   'gui_OutputFcn',  @mbc_vmfmm_OutputFcn, ...
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


% --- Executes just before mbc_vmfmm is made visible.
function mbc_vmfmm_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to mbc_vmfmm (see VARARGIN)

% Choose default command line output for mbc_vmfmm
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes mbc_vmfmm wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = mbc_vmfmm_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in load_data_GT.
function load_data_GT_Callback(hObject, eventdata, handles)
% hObject    handle to load_data_GT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Select the appropriate file
[FileName,PathName,FilterIndex] = uigetfile;

% Load samples and labels
load(strcat(PathName, FileName));

% global data 
handles.vmfSample = vmfSample;
handles.labels = labels;
guidata(hObject,handles);

% Display data in sphere with label
spread_gui(handles.axes1, vmfSample', labels);

% Enable clustering button
set(handles.cluster_data,'Enable','on');

% --- Executes on key press with focus on load_data_GT and none of its controls.
function load_data_GT_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to load_data_GT (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over load_data_GT.
function load_data_GT_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to load_data_GT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in load_data_only.
function load_data_only_Callback(hObject, eventdata, handles)
% hObject    handle to load_data_only (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Select the appropriate file
[FileName,PathName,FilterIndex] = uigetfile;

% Load samples
load(strcat(PathName, FileName));

% global data 
handles.vmfSample = vmfSample;
handles.labels = [];
guidata(hObject,handles);

% Display data in sphere without label
spread_gui(handles.axes1, vmfSample');

% Enable clustering button
set(handles.cluster_data,'Enable','on');

% --- Executes on button press in load_data_GT.
function load_data_GT_Callback_Callback(hObject, eventdata, handles)
% hObject    handle to load_data_GT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in load_data_only.
function load_data_only_Callback_Callback(hObject, eventdata, handles)
% hObject    handle to load_data_only (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in save_plot.
function save_plot_Callback(hObject, eventdata, handles)
% hObject    handle to save_plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Select file name
FileName = uiputfile('*.*');
export_fig(handles.axes1, strcat(FileName, '.png'));


% --- Executes during object creation, after setting all properties.
function figure1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
addpath('exportfig');
addpath('vmfmatlab');


% --- Executes on button press in gen_vmfmm.
function gen_vmfmm_Callback(hObject, eventdata, handles)
% hObject    handle to gen_vmfmm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Select the appropriate file
[FileName,PathName,FilterIndex] = uigetfile;

% Load the mixture model
load(strcat(PathName, FileName));

% Generate samples from the mixture model
numSamp = 10000;
[vmfSample,labels] = emsamp(mixture, numSamp);

% global data 
handles.vmfSample = vmfSample;
handles.labels = labels;
guidata(hObject,handles);

% Display data in sphere with label
spread_gui(handles.axes1, vmfSample', labels);

% Enable clustering button
set(handles.cluster_data,'Enable','on');

% --- Executes on button press in cluster_data.
function cluster_data_Callback(hObject, eventdata, handles)
% hObject    handle to cluster_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Load necessary data from global data 
handles = guidata(hObject);
vmfSample = handles.vmfSample;

% Perform clustering data
kMax = 15;

% get model parameters with k_max
params = bd_vmfmm(vmfSample, kMax);
params = annihilateComp(params); % eliminate empty clusters

% get model parameters with k_max
[~, allIC, allClust] = get_HC_IC_BD(vmfSample, params); 

% Model selection
y = allIC.BIC;
x = 1: length(y);

% Select the number of clusters automatically (s.t. k_max number of
% clusters)
[numComp,~, ~] = wplr(x,y,[1 30]); % WPLR-30

% Get the final clustering results
finalClust = allClust(:, numComp);

% Display data in sphere with label
spread_gui(handles.axes1, vmfSample', finalClust);

% Disable clustering button
set(handles.cluster_data,'Enable','off');
