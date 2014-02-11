function varargout = ImageBlurrer(varargin)
% IMAGEBLURRER MATLAB code for ImageBlurrer.fig
%      IMAGEBLURRER, by itself, creates a new IMAGEBLURRER or raises the existing
%      singleton*.
%
%      H = IMAGEBLURRER returns the handle to a new IMAGEBLURRER or the handle to
%      the existing singleton*.
%
%      IMAGEBLURRER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in IMAGEBLURRER.M with the given input arguments.
%
%      IMAGEBLURRER('Property','Value',...) creates a new IMAGEBLURRER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ImageBlurrer_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ImageBlurrer_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ImageBlurrer

% Last Modified by GUIDE v2.5 19-Apr-2012 17:53:41

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ImageBlurrer_OpeningFcn, ...
                   'gui_OutputFcn',  @ImageBlurrer_OutputFcn, ...
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


% --- Executes just before ImageBlurrer is made visible.
function ImageBlurrer_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ImageBlurrer (see VARARGIN)

% Choose default command line output for ImageBlurrer
handles.output = hObject;

% Determine the size of the MatLab parallel processing pool
handles.parallel_CPUs = matlabpool('size');
if handles.parallel_CPUs > 0
    set(handles.parallel_check,'Value',1)    
    set(handles.CPU_edit,'String',num2str(handles.parallel_CPUs))
else
    set(handles.parallel_check,'Value',0)
    set(handles.CPU_edit,'String','0')
end

% Initial path from which source and target files are searched
handles.source_path = pwd;
handles.target_path = pwd;


% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ImageBlurrer wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ImageBlurrer_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pick_source_button.
function pick_source_button_Callback(hObject, eventdata, handles)
% hObject    handle to pick_source_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Request folder containing video files to be blurred
handles.source_path = uigetdir(handles.source_path, ...
    'Pick source root folder');

set(handles.source_edit,'String',handles.source_path);

handles.source_full_path = get(handles.source_edit,'String');

guidata(hObject,handles)


% --- Executes on button press in pick_target_button.
function pick_target_button_Callback(hObject, eventdata, handles)
% hObject    handle to pick_target_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Request target folder to save analysis results in
handles.target_path = uigetdir(handles.target_path, ...
    'Pick target root folder');

set(handles.target_edit,'String',handles.target_path);

guidata(hObject,handles)


function source_edit_Callback(hObject, eventdata, handles)
% hObject    handle to source_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of source_edit as text
%        str2double(get(hObject,'String')) returns contents of source_edit as a double

set(handles.source_edit,'String',handles.source_full_path)

guidata(hObject,handles)



% --- Executes during object creation, after setting all properties.
function source_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to source_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function target_edit_Callback(hObject, eventdata, handles)
% hObject    handle to target_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of target_edit as text
%        str2double(get(hObject,'String')) returns contents of target_edit as a double

set(handles.target_edit,'String',handles.target_path)

guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function target_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to target_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in parallel_check.
function parallel_check_Callback(hObject, eventdata, handles)
% hObject    handle to parallel_check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of parallel_check

handles.parallel_CPUs = matlabpool('size');

if handles.parallel_CPUs > 0
    matlabpool close
else
    if str2double(get(handles.CPU_edit,'String')) > 0
        matlabpool(str2double(get(handles.CPU_edit,'String')))
    else
        matlabpool
    end
end
    
handles.parallel_CPUs = matlabpool('size');
if handles.parallel_CPUs > 0
    set(handles.parallel_check,'Value',1)
else
    set(handles.parallel_check,'Value',0)
end
set(handles.CPU_edit,'String',num2str(handles.parallel_CPUs))


function CPU_edit_Callback(hObject, eventdata, handles)
% hObject    handle to CPU_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CPU_edit as text
%        str2double(get(hObject,'String')) returns contents of CPU_edit as a double

handles.CPUs = ceil(str2double(get(handles.CPU_edit,'String')));
set(handles.CPU_edit,'String',num2str(handles.CPUs))

% --- Executes during object creation, after setting all properties.
function CPU_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CPU_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in start_button.
function start_button_Callback(hObject, eventdata, handles)
% hObject    handle to start_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

source_exist_value = ...
    exist(get(handles.source_edit,'String')','file');
target_exist_value = ...
    exist(get(handles.target_edit,'String')','file');

if ~source_exist_value
    msgbox('Choose valid source path.','Invalid source path','error')
    return
elseif ~target_exist_value
    msgbox('Choose valid target path','Invalid target path','error')
    return
end

% ------------
% Create input cell

target_root_directory = ...
    get(handles.target_edit,'String');
% --
% Create input batch that will be used in EvaluateVideo function call

%Get the source and target root directory
source_root_directory = ...
    get(handles.source_edit,'String');
%Call function that finds the paths to the
[avi_file_paths,target_paths] = ...
    findvideos(source_root_directory,target_root_directory);

% Batch to supply to EvaluateVideo.m function
videos = {avi_file_paths{:};target_paths{:}}.';


% ------------
% Create structure containing the analysis parameters
parameters = struct;

%Micrometers per pixel in the video
parameters.filtersize = ...
    str2double(get(handles.filtersize_edit,'String'));

%Raw frames merged into one compressed frame
parameters.filtersigma = ...
    str2double(get(handles.filtersigma_edit,'String'));

% -------------------
%Call to blurring function

fprintf('Call to blurring function...\n')
block_box = msgbox(['Do not use the interface while the blurring is running.' ...
    ' You might accidentally start a second run.' ...
    ' Check command line for progress.'],...
    'Analysis started',...
    'warn');

tic

BlurVideo(videos,parameters);

try
    close(block_box)
catch ME
    % bla
end
elapsed_time = toc;
fprintf('Blurring completed, computation time %f seconds\n',...
    elapsed_time)

guidata(hObject,handles)



% --- Executes on button press in save_parameters_button.
function save_parameters_button_Callback(hObject, eventdata, handles)
% hObject    handle to save_parameters_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% ------------
% Create structure containing the analysis parameters
parameters = struct;

%Size of blur filter
parameters.filtersize = ...
    str2double(get(handles.filtersize_edit,'String'));

%Sharpness of blur filter
parameters.filtersigma = ...
    str2double(get(handles.filtersigma_edit,'String'));

[save_file,save_path] = uiputfile(pwd,'Save parameters in which file?');
save([save_path filesep save_file],'parameters')

guidata(hObject,handles)


% --- Executes on button press in load_parameter_button.
function load_parameter_button_Callback(hObject, eventdata, handles)
% hObject    handle to load_parameter_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% ------------
% Load parameters from file and assign to edit fields in the GUI

% ----
% Load parameters from file

%Which parameter file should be loaded?
[load_file,load_path] = uigetfile(pwd,'Select parameter file');

% Load from chosen file
load([load_path filesep load_file])

% ----
%Assign values from parameter file to edit fields in the GUI

%Micrometers per pixel in the video
set(handles.filtersize_edit,'String', ...
    num2str(parameters.filtersize))

%Raw frames merged into one compressed frame
set(handles.filtersigma_edit,'String', ...
    num2str(parameters.filtersigma))

guidata(hObject,handles)





function filtersize_edit_Callback(hObject, eventdata, handles)
% hObject    handle to filtersize_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


%Make sure that absolute length change analysis parameter is not negative
filtersize_value = ...
    str2int(get(handles.abs_length_change_edit,'String'));

% Make zero when negative, leave at same value if positive
filtersize_value = ...
    (round(filtersize_value)>=3).*round(filtersize_value) ...
    + 3.*(round(filtersize_value)<3);

%Reassign to edit field
set(handles.filtersize_edit,'String', ...
    num2str(filtersize_value))


% Hints: get(hObject,'String') returns contents of filtersize_edit as text
%        str2double(get(hObject,'String')) returns contents of filtersize_edit as a double

guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function filtersize_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to filtersize_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function filtersigma_edit_Callback(hObject, eventdata, handles)
% hObject    handle to filtersigma_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%Make sure that absolute length change analysis parameter is not negative

filtersigma_value = ...
    str2int(get(handles.filtersigma_edit,'String'));

% Make zero when negative, leave at same value if positive
filtersigma_value = ...
    (round(filtersigma_value)>0).*round(filtersigma_value) ...
    + 0.25.*(round(filtersigma_value)<=0);

%Reassign to edit field
set(handles.filtersigma_edit,'String', ...
    num2str(filtersigma_value))

% Hints: get(hObject,'String') returns contents of filtersize_edit as text
%        str2double(get(hObject,'String')) returns contents of filtersize_edit as a double

guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function filtersigma_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to filtersigma_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
