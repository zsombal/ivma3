function varargout = TraceDropper(varargin)
% TRACEDROPPER MATLAB code for TraceDropper.fig
%      TRACEDROPPER, by itself, creates a new TRACEDROPPER or raises the existing
%      singleton*.
%
%      H = TRACEDROPPER returns the handle to a new TRACEDROPPER or the handle to
%      the existing singleton*.
%
%      TRACEDROPPER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TRACEDROPPER.M with the given input arguments.
%
%      TRACEDROPPER('Property','Value',...) creates a new TRACEDROPPER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before TraceDropper_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to TraceDropper_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help TraceDropper

% Last Modified by GUIDE v2.5 16-May-2014 16:07:36

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @TraceDropper_OpeningFcn, ...
                   'gui_OutputFcn',  @TraceDropper_OutputFcn, ...
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


% --- Executes just before TraceDropper is made visible.
function TraceDropper_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to TraceDropper (see VARARGIN)

% Choose default command line output for TraceDropper
handles.output = hObject;

% ----
% Set properties of the GUI elements at time of initialization
% No ticks for trace display axis
set(handles.trace_display,'XTick',[],'YTick',[])
set(handles.trace_display,'Box','on')
% No ticks for progress display axis
set(handles.progress_display,'XTick',[],'YTick',[],...
    'Box','on')
% No ticks for filament video display axis
set(handles.frame_video_display,'XTick',[],'YTick',[],...
    'Box','on')

% Make the text displaying the batch file number invisible
set(handles.batch_number_text,'Visible','off')

% Set the write_plain_text checkbutton variable and also set the checkbox
% to unchecked
handles.write_plain_text = 1;
set(handles.write_text_checkbox,'Value',1)

% Set colormap for image plots
colormap('gray')

% Put keyboard commands into the GUI
set(handles.keyboard_commands_text,'String', ...
    sprintf(['Keyboard commands:\n' ...
    'Return - next trace\n' ...
    'Backspace - drop and hide\n', ...
    'Shift - keep and show\n', ...
    'Space - keep and hide\n', ...
    'Control - drop and show\n']))

% Instantiate cell array to keep track of trace choices - for UNDO function
handles.undo_cell = {};

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes TraceDropper wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = TraceDropper_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in choose_folder.
function choose_folder_Callback(hObject, eventdata, handles)
% hObject    handle to choose_folder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isfield(handles,'current_path')
    handles.current_path = uigetdir(handles.current_path);
else
    handles.current_path = uigetdir;
end

% Set batch progress to 0 to indicate that we are working in single result
% folder mode, not in batch mode
handles.batch_progress = 0;

handles = load_result_folder(hObject,handles);

% Update handles structure
guidata(hObject, handles);



% --- Executes on button press in choose_batch.
function choose_batch_Callback(hObject, eventdata, handles)
% hObject    handle to choose_batch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[batch_file,batch_path] = ...
    uigetfile('*.mat','Pick a TraceDropper_listfile.');
batch_saveload = load([batch_path filesep batch_file]);
handles.batch = batch_saveload.target_folder_list;
handles.batch_progress = 0;

% Flag indicating if batch analysis has been started set to false
handles.batch_started = false;

handles = next_folder_in_batch(hObject,handles);

% Update handles structure
guidata(hObject, handles);


function current_folder_edit_Callback(hObject, eventdata, handles)
% hObject    handle to current_folder_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of current_folder_edit as text
%        str2double(get(hObject,'String')) returns contents of current_folder_edit as a double

% Make sure that user entry does not change path, so restore after user has
% entered different letters
if isfield(handles,'current_path')
    set(handles.current_folder_edit,'String',...
        handles.current_path)
else
    set(handles.current_folder_edit,'String',...
        'Current folder will appear here')
end

% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function current_folder_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to current_folder_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on key press with focus on figure1 and none of its controls.
function figure1_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)

key_string = eventdata.Key;
if ismember(key_string,{'return','backspace','space','control','shift'})
    if strcmp(key_string,'shift')
        handles = put_this_trace(hObject,handles,'keepshow');
    elseif strcmp(key_string,'backspace')
        handles = put_this_trace(hObject,handles,'drophide');
    elseif strcmp(key_string,'space')
        handles = put_this_trace(hObject,handles,'keephide');
    elseif strcmp(key_string,'control')
        handles = put_this_trace(hObject,handles,'dropshow');
    elseif strcmp(key_string,'return')
        handles = put_this_trace(hObject,handles,'next');
    end
end

guidata(hObject,handles)



% --- Executes on button press in drop_all_button.
function drop_all_button_Callback(hObject, eventdata, handles)
% hObject    handle to drop_all_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


handles.keep_vector(:) = 0;
handles = update_plots(hObject,handles);
handles.undo_cell = {};

guidata(hObject,handles)



% --- Executes on button press in keep_all_button.
function keep_all_button_Callback(hObject, eventdata, handles)
% hObject    handle to keep_all_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.keep_vector(:) = 1;
handles = update_plots(hObject,handles);
handles.undo_cell = {};

guidata(hObject,handles)



% --- Executes on button press in show_all_button.
function show_all_button_Callback(hObject, eventdata, handles)
% hObject    handle to show_all_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.show_vector(:) = 1;
handles = update_plots(hObject,handles);
handles.undo_cell = {};

guidata(hObject,handles)


% --- Executes on mouse press over axes background.
function VL_diagram_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to VL_diagram (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Pick up current mouse pointer position from the (V,L) display
current_point = get(handles.VL_diagram,'CurrentPoint');
% Extract horizontal coordinate and make sure it is positive
handles.clicked_length = (current_point(1,1));
if handles.clicked_length<0
    handles.clicked_length = 0;
end
% Extract the vertical coordinate and make sure it is positive
handles.clicked_velocity = (current_point(1,2));
if handles.clicked_velocity<0
    handles.clicked_velocity = 0;
end
%fprintf('Clicked velocity: %f\n',handles.clicked_velocity)
%fprintf('Clicked length: %f\n',handles.clicked_length)


% Calculate the distance to all (V,L) points
distance_vector = ...
    (abs((handles.clicked_length...
    -handles.current_average_lengths(handles.show_vector))./ ...
    handles.current_max_average_length) ...
    + abs((handles.clicked_velocity...
    -handles.current_trace_velocities(handles.show_vector))./ ...
    handles.current_max_trace_velocity)).';

% Find (V,L) point with minimum distance to put the view on the trace
% corresponding to the (V,L) that is the closest to the clicked point
[~,min_index] = min(distance_vector);
show_indices = find(handles.show_vector);
handles.trace_in_view = show_indices(min_index);

handles = update_plots(hObject,handles);

guidata(hObject,handles)


% --- Executes on mouse press over axes background.
function trace_display_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to trace_display (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Pick up current mouse pointer position from the trace display
current_point = get(handles.trace_display,'CurrentPoint');
% Extract horizontal coordinate and make sure it is inside the video size
% range
handles.trace_xx = round(current_point(1,1));
if handles.trace_xx<0
    handles.trace_xx = 0;
elseif handles.trace_xx>handles.video_width
    handles.trace_xx = handles.video_width;
end
% Extract the vertical coordinate and make sure it is inside the video size
% range
handles.trace_yy = round(current_point(1,2));
if handles.trace_yy<0
    handles.trace_yy = 0;
elseif handles.trace_yy>handles.video_height
    handles.trace_yy = handles.video_height;
end

% Bring the clicked trace into view
clicked_index = handles.index_matrix(handles.trace_yy,handles.trace_xx);
if clicked_index > 0
    handles.trace_in_view = clicked_index;
end

%Update the GUI display accordingly
handles = update_plots(hObject,handles);

guidata(hObject,handles)


% --- Executes on button press in save_button.
function save_button_Callback(hObject, eventdata, handles)
% hObject    handle to save_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = save_keepdrop_mask(hObject,handles);

if isfield(handles,'batch') && handles.batch_progress>0
    % Call function that executes the jump to the next result folder in the
    % batch
    handles = next_folder_in_batch(hObject,handles);
end

% Update handles structure
guidata(hObject,handles)


% --- Executes on button press in write_text_checkbox.
function write_text_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to write_text_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of write_text_checkbox

handles.write_plain_text = get(handles.write_text_checkbox,'Value');
guidata(hObject,handles)


% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = undo_last(hObject,handles);

guidata(hObject,handles)


% --- Can be called from other functions
function handles = load_result_folder(hObject,handles)
% handles    empty - handles not created until after all CreateFcns called
% hObject    handle to object that called this function
%
% This function loads the results from the currently chosen result folder


% Check if the result folder to be loaded does exist on disk
if ~exist([handles.current_path filesep 'Analysis_Results.mat'],'file')
    % If result folder does not exist on disk
    handles.batch_started = false;
    fprintf('Non-existent result folder skipped.\n')
    % Indicate that this directory should not be treated by the
    % TraceDropper, here for the reason that it simply does not exist
else
    % If result folder to be load does exist on disk, load it into the
    % TraceDropper for the user to remove bad traces
    
    % Display current path in the current_folder_edit box
    set(handles.current_folder_edit,'String',handles.current_path)
    
    % If in batch mode, display current file number
    if handles.batch_progress > 0
        set(handles.batch_number_text,'Visible','on',...
            'String',...
            sprintf('Batch file %d/%d',handles.batch_progress,...
            numel(handles.batch)))
    else
        set(handles.batch_number_text,'Visible','off')
    end
    
    
    % Load the current folder's results and initialize containers and
    % bookkeeping arrays in the handles object
    handles.current_results = ...
        load([handles.current_path filesep 'Analysis_Results.mat']);
    handles.video_width = ...
        handles.current_results.trace_results.video_properties.width;
    handles.video_height = ...
        handles.current_results.trace_results.video_properties.height;
    handles.current_tracenumber = ...
        handles.current_results.trace_results.video_properties.traces;
    %For drawing of traces at different grey values
    rand_column = 0.5 + 0.3.*rand(handles.current_tracenumber,1);
    handles.random_colormap = cat(1,...
        cat(2,rand_column,rand_column,rand_column), ...
        [1,0,0],[1 1 1],[0 0 0]);
    
    % contains true for the traces that have not been dropped so far
    handles.keep_vector = ...
        true(1,handles.current_tracenumber);
    % Contains true for the traces that are still to be shown in the plots
    handles.show_vector = ...
        true(1,handles.current_tracenumber);
    
    %Load the linear indices of pixels contained in the individual traces
    handles.trace_indices = ...
        {handles.current_results.trace_results.trace_results. ...
        trace_pixel_indices_linear};
    handles.individual_filament_indices = { ...
        handles.current_results.trace_results.trace_results. ...
        filament_pixel_indices_linear};
    
    % Load the frame-to-frame velocity data from the current result folder
    handles.frame_to_frame_velocities = { ...
        handles.current_results.trace_results.trace_results.frame_to_frame_velocities};
    
    % Load the (V,L) data from the current result folder
    handles.current_average_lengths = ...
        real([handles.current_results.trace_results.trace_results.average_filament_length]);
    handles.current_max_average_length = ...
        max(handles.current_average_lengths);
    handles.current_trace_velocities = ...
        real([handles.current_results.trace_results.trace_results.trace_velocity]);
    handles.current_max_trace_velocity = ...
        max(handles.current_trace_velocities);
    % Load the length time traces from the current result folder
    handles.current_length_timetraces = ...
        {handles.current_results.trace_results.trace_results.filament_lengths};
    handles.current_length_timetraces = cellfun( ...
        @(xx) real(xx),handles.current_length_timetraces,'UniformOutput',false);
    % Load the solidity time traces from the current result folder
    handles.current_solidity_timetraces = ...
        {handles.current_results.trace_results.trace_results.filament_solidities};
    
    %Load the bounding boxes for filament traces
    handles.current_bounding_box = ...
        {handles.current_results.trace_results.trace_results.trace_bounding_box};
    
    % Determine delta t between two frames and the present frame number for
    % each of the tracked filaments in each of the traces
    handles.current_delta_t = ...
        handles.current_results.trace_results.video_properties.frame_rate.^(-1);
    handles.current_present_frames = ...
        [handles.current_results.trace_results.trace_results.present_frames];
    
    % Has the current result file been treated with the TraceDropper and
    % therefore holds a tracedrop_mask? If not, it will have to be treated now,
    % so set the flag batch_started to yes.
    if isfield(handles.current_results,'tracedrop_mask')
        handles.batch_started = false;
    else
        handles.batch_started = true;
    end
    
    % Initialize the trace view on the trace with index one
    handles.trace_in_view = 1;
    
    handles = update_plots(hObject,handles);
    
    % Clear the UNDO history vector
    handles.undo_cell = {};
    
end

guidata(hObject,handles)



% --- Can be called from other functions
function handles = update_plots(hObject,handles)
% handles    empty - handles not created until after all CreateFcns called
% hObject    handle to object that called this function
%
% This function redraws the displays according to current status of the
% filament dropper


% ----
% Plotting of the frame-to-frame velocity histogram

hist(handles.f2f_histogram,[handles.frame_to_frame_velocities{ ...
    handles.keep_vector}],25)
xlabel(handles.f2f_histogram,'V_{f2f}[\mum/s]')
ylabel(handles.f2f_histogram,'Count')

% ----
% Plotting of the trace display

%Matrix to contain all traces' indices to the random color map
handles.colorindex_matrix = ...
    (handles.current_tracenumber+2) ...
    .*ones(handles.video_height,handles.video_width,'uint64');
%Matrix to contain the indices of all traces for point-and-click selection
%of filaments
handles.index_matrix = ...
    zeros(handles.video_height,handles.video_width,'uint64');


% Add all traces that are in show vector into the grey matrix for being
% displayed. Also add to the index matrix for later point-and-click
for tt = find(handles.show_vector)    
    handles.colorindex_matrix(handles.trace_indices{tt}) = tt;
    handles.index_matrix(handles.trace_indices{tt}) = tt;
end


% ----
% Plotting of the (V,L) diagram
cla(handles.VL_diagram)
quick_handle = plot(handles.VL_diagram,...
    handles.current_average_lengths(handles.show_vector),...
    handles.current_trace_velocities(handles.show_vector),...
    'ko','MarkerSize',6);
set(quick_handle,'HitTest','off')
set(handles.VL_diagram,'NextPlot','add');
xlabel(handles.VL_diagram,'Filament length[\mum]')
ylabel(handles.VL_diagram,'Filament velocity[\mum/s]')
set(handles.VL_diagram,...
    'XLim',[0 1.05.*max(handles.current_average_lengths)],...
    'YLim',[0 1.05.*max(handles.current_trace_velocities)])

% ---- Things that can be plotted when a specific filament is currently in
% view

if handles.trace_in_view
    
    % Add the trace currently in view as red trace and let the other
    % traces disappear into the background of the trace display
    handles.colorindex_matrix( ...
        handles.trace_indices{handles.trace_in_view}) = ...
        handles.current_tracenumber+1;
    handles.index_matrix(handles.trace_indices{handles.trace_in_view}) = ...
        handles.trace_in_view;
    
    % Draw a blow-up of the current trace in the frame_video_display axes
    this_trace_image_matrix = ...
        double(handles.current_results.trace_results. ...
        trace_results(handles.trace_in_view).trace_image);
    this_trace_image_matrix(this_trace_image_matrix==1) = ...
        handles.current_tracenumber + 3;
    this_trace_image_matrix(this_trace_image_matrix==0) = ...
        handles.current_tracenumber + 2;
    image(this_trace_image_matrix,'Parent',handles.frame_video_display)
    set(handles.frame_video_display,...
        'XTick',[],'YTick',[])
    
    
    
    
    % Draw a crosshair for the (V,L) point of the filament that is currently in
    % view    
    quick_handle = plot(handles.VL_diagram,...
        [0 1.05.*max(handles.current_average_lengths)],...
        handles.current_trace_velocities(handles.trace_in_view).*ones(1,2),...
        'k-');
    set(quick_handle,'HitTest','off')
    quick_handle = plot(handles.VL_diagram,...
        handles.current_average_lengths(handles.trace_in_view).*ones(1,2),...
        [0 1.05.*max(handles.current_trace_velocities)],...
        'k-');
    set(quick_handle,'HitTest','off')
    % Draw a filled square at the (V,L) data point for the filament
    % currently in view
    quick_handle = plot(handles.VL_diagram,...
        handles.current_average_lengths(handles.trace_in_view),...
        handles.current_trace_velocities(handles.trace_in_view),...
        'ks','MarkerSize',8,'MarkerFaceColor',[0 0 0]);
    set(quick_handle,'HitTest','off')
    set(handles.VL_diagram,'NextPlot','replace');
    
%     % ----
%     % Plotting of the length and solidity time traces
%     axes(handles.time_course_display);
%     time_support = ...
%         ([1:handles.current_present_frames(handles.trace_in_view)]-0.5) ...
%         .* handles.current_delta_t;
%     [axis_handles,line_handle_one,line_handle_two] = ...
%         plotyy( ...
%         time_support,...
%         handles.current_length_timetraces{handles.trace_in_view}, ...
%         time_support, ...
%         handles.current_solidity_timetraces{handles.trace_in_view});
%     ylabel(axis_handles(1),'Fil. length[\mum]')
%     xlabel(axis_handles(1),'Time[s]')
%     ylabel(axis_handles(2),'Solidity')
%     
%     set(axis_handles(1),'YColor',[0 0 0])
%     set(axis_handles(1),'XTick',[time_support(1) time_support(end)],...
%         'XLim',[time_support(1),time_support(end)],...
%         'YLim',[0 ...
%         1.05.*max(handles.current_length_timetraces{handles.trace_in_view})],...
%         'YTick',linspace(0, ...
%         1.05.*max(handles.current_length_timetraces{handles.trace_in_view}),...
%         5),...
%         'TickLength',[0 0])
%     set(axis_handles(2),'XTick',[],...
%         'XLim',[time_support(1),time_support(end)],...
%         'YLim',[0 1],...
%         'YTick',linspace(0,1,5),...
%         'TickLength',[0 0])
%     set(axis_handles(2),'YColor',[1 0 0])
%     set(line_handle_one,'Color',[0 0 0])
%     set(line_handle_two,'Color',[1 0 0],...
%         'LineStyle','--')
    
    % ----
    % Plotting the sorting progress and sorting classes
    plot(handles.progress_display,[0.5 0.5],[0 handles.trace_in_view],...
        'k-','LineWidth',8)
    set(handles.progress_display,'NextPlot','add')
    plot(handles.progress_display,...
        handles.keep_vector,[1:handles.current_tracenumber],...
        'ko')
    no_show_indices = find(~handles.show_vector);
    if ~isempty(no_show_indices)
        plot(handles.progress_display,...
            [-0.5 1.5],cat(1,no_show_indices,no_show_indices),...
            'k-')
    end
    set(handles.progress_display,'YLim',[0 handles.current_tracenumber])
    set(handles.progress_display,'XLim',[-0.5 1.5])
    set(handles.progress_display,'YTick',[],'XTick',[0 1],...
        'TickLength',[0 0])
    set(handles.progress_display,'NextPlot','replace')

end

colormap(handles.random_colormap);
handles.trace_image = ...
    image(handles.colorindex_matrix,'Parent',handles.trace_display);
% Add a bounding box rectangle for the trace that is in view
if handles.trace_in_view
    set(handles.trace_display,'NextPlot','add')
    rectangle( ...
        'Position',handles.current_bounding_box{handles.trace_in_view}, ...
        'Parent',handles.trace_display, ...
        'HitTest','off', ...
        'LineWidth',2,'EdgeColor',[0 0 0])
    set(handles.trace_display,'NextPlot','replace')
end

% Assign button down functions for the axes that can be clicked for trace
% selection
set(handles.VL_diagram,'ButtonDownFcn',...
    @(hObject,eventdata) ...
    TraceDropper('VL_diagram_ButtonDownFcn', ...
    hObject, eventdata, guidata(hObject)));
set(handles.trace_image,'HitTest',...
    'off')
set(handles.trace_display,'ButtonDownFcn', ...
    @(hObject,eventdata) ...
    TraceDropper('trace_display_ButtonDownFcn', ...
    hObject, eventdata, guidata(hObject)))

set(handles.trace_display,'XTick',[],'YTick',[])

% ---- Can be called by callback functions
function handles = put_this_trace(hObject,handles,put_where)
% handles    empty - handles not created until after all CreateFcns called
% put_where put_where - string that determines into which class the
%           current trace should be put
%
% This function sorts the current trace into the class indicated by
% put_where, and adjust the show_vector (true for traces that are still to
% be displayed in the trace dropper interface) and the keep_vector (true
% for traces that are kept in the result data set for further analysis)
% accordingly

last_state = [handles.trace_in_view, ...
    handles.keep_vector(handles.trace_in_view), ...
    handles.show_vector(handles.trace_in_view)];

handles.undo_cell = [handles.undo_cell,last_state];

if strcmp(put_where,'drophide')
    % Drop the current trace from results and remove from displays
    handles.keep_vector(handles.trace_in_view) = 0;
    handles.show_vector(handles.trace_in_view) = 0;
elseif strcmp(put_where,'keephide')
    % Keep the current trace, but remove from the displays
    handles.keep_vector(handles.trace_in_view) = 1;
    handles.show_vector(handles.trace_in_view) = 0;
elseif strcmp(put_where,'dropshow')
    % Drop the current trace, but keep for displays
    handles.keep_vector(handles.trace_in_view) = 0;
    handles.show_vector(handles.trace_in_view) = 1;
elseif strcmp(put_where,'keepshow')
    % Keep the current trace, and keep for displays
    handles.keep_vector(handles.trace_in_view) = 1;
    handles.show_vector(handles.trace_in_view) = 1;
end

if sum(handles.show_vector) > 0
    
    
    % If there are traces left that could potentially be dropped, keep
    % going through the filaments that are still visible
    if handles.trace_in_view == handles.current_tracenumber
        handles.trace_in_view = 1;
    else
        handles.trace_in_view = handles.trace_in_view + 1;
    end
    while ~handles.show_vector(handles.trace_in_view)
        if handles.trace_in_view == handles.current_tracenumber
            handles.trace_in_view = 1;
        else
            handles.trace_in_view = handles.trace_in_view + 1;
        end
    end
    
elseif isfield(handles,'batch') && handles.batch_progress>0
    % Call function that executes the jump to the next result folder in the
    % batch
    handles = next_folder_in_batch(hObject,handles);
end

handles = update_plots(hObject,handles);

guidata(hObject,handles)


% ---- Can be called by callback functions
function handles = undo_last(hObject,handles)
% handles    empty - handles not created until after all CreateFcns called
% put_where put_where - string that determines into which class the
%           current trace should be put
%
% This function restores the last saved values in the undo cell

% Retrieve and remove last action

last_action = false;

if numel(handles.undo_cell)>0
        
    last_action = handles.undo_cell{end};
    
    % Remove last action from undo cell
    if numel(handles.undo_cell)==1
        handles.undo_cell = {};
    elseif numel(handles.undo_cell)>1
        handles.undo_cell = handles.undo_cell(1:end-1);
    end
end

% Restore filament to before last action

if last_action
   
    handles.keep_vector(last_action(1)) = last_action(2);
    handles.show_vector(last_action(1)) = last_action(3);
    handles.trace_in_view = last_action(1);
    
end

handles = update_plots(hObject,handles);


% ---- Can be called from other functions
function handles = next_folder_in_batch(hObject,handles)
% hObject    handle to trace_display (see GCBO)
% handles    structure with handles and user data (see GUIDATA)
%
% Sets the next folder in the folder batch as the current path, stores the
% results from the folder worked on before with the existen results in that
% folder

% Store trace dropping results
if handles.batch_progress>0 && handles.batch_started
    handles = save_keepdrop_mask(hObject,handles);
end

% Return message when the batch is finished
if handles.batch_progress == numel(handles.batch)
    msgbox('Batch finished, load next batch or result folder to continue.')
end

% Jump to next results folder in batch that has not yet been treated with
% the TraceDropper
if handles.batch_progress < numel(handles.batch)
    handles.batch_progress = handles.batch_progress + 1;
    handles.current_path = handles.batch{handles.batch_progress};
    handles = load_result_folder(hObject,handles);
end
while handles.batch_progress < numel(handles.batch)  && ...
        ~handles.batch_started
    handles.batch_progress = handles.batch_progress + 1;
    handles.current_path = handles.batch{handles.batch_progress};
    handles = load_result_folder(hObject,handles);
end


% ---- Save the current keepdrop mask with the existent analysis results
function handles = save_keepdrop_mask(hObject,handles)
% hObject    handle to trace_display (see GCBO)
% handles    structure with handles and user data (see GUIDATA)
%
% Saves the current keepdrop vector with the existen analysis results

% Append tracedrop_mask to existent analysis results
tracedrop_mask = handles.keep_vector;
save([handles.current_path filesep 'Analysis_Results.mat'], ...
    'tracedrop_mask','-append');


% ----
% Write to plain text file if requested by options

if handles.write_plain_text
    
    % Trace properties
    
    % Tab delimited file with all trace properties
    trace_text_file = ...
        [handles.current_path filesep 'Dropped_Trace_Results.dat'];
    names_of_fields = ...
        {'average_filament_length',...
        'average_filament_width',...
        'trace_velocity',...
        'trace_solidity'};
    field_headers = ...
        {'Length[mum]',...
        'Width[mum]',...
        'Velocity[mum/s]',...
        'Solidity'};
    
    number_of_columns = numel(names_of_fields);
    
    % Write the headers for data colums
    trace_result_file_writeobject = ...
        fopen(trace_text_file,'w');
    for ff = 1:number_of_columns
        fprintf(trace_result_file_writeobject,...
            '%s\t',field_headers{ff});
    end
    fprintf(trace_result_file_writeobject,'\n');
    fclose(trace_result_file_writeobject);
    
    % Compile write array and write to file
    
    kept_traces = find(handles.keep_vector);
    if ~isempty(kept_traces)
        write_array = zeros(numel(kept_traces),number_of_columns);
        for ff = 1:number_of_columns
            write_array(:,ff) = ...
                real([handles.current_results.trace_results. ...
                trace_results(kept_traces).(names_of_fields{ff})]);
        end
        dlmwrite(trace_text_file,write_array,'-append','delimiter','\t')
    end
    
end


