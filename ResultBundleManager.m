function varargout = ResultBundleManager(varargin)
% RESULTBUNDLEMANAGER MATLAB code for ResultBundleManager.fig
%      RESULTBUNDLEMANAGER, by itself, creates a new RESULTBUNDLEMANAGER or raises the existing
%      singleton*.
%
%      H = RESULTBUNDLEMANAGER returns the handle to a new RESULTBUNDLEMANAGER or the handle to
%      the existing singleton*.
%
%      RESULTBUNDLEMANAGER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in RESULTBUNDLEMANAGER.M with the given input arguments.
%
%      RESULTBUNDLEMANAGER('Property','Value',...) creates a new RESULTBUNDLEMANAGER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ResultBundleManager_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ResultBundleManager_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ResultBundleManager

% Last Modified by GUIDE v2.5 24-Jan-2012 15:34:17

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ResultBundleManager_OpeningFcn, ...
                   'gui_OutputFcn',  @ResultBundleManager_OutputFcn, ...
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


% --- Executes just before ResultBundleManager is made visible.
function ResultBundleManager_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ResultBundleManager (see VARARGIN)

% Choose default command line output for ResultBundleManager
handles.output = hObject;

% Initialize variables to keep track of where the user is working
handles.bundle_path = false;
handles.source_path = false;

% Initialize variables to keep track of user entries to the keyword edit
% fields
handles.any_keywords_edited = false;
handles.any_keywords = {};
handles.with_keywords_edited = false;
handles.with_keywords = {};
handles.not_keywords_edited = false;
handles.not_keywords = {};


% Initialize variable to contain histogram binning V_high
handles.V_high = str2double(get(handles.V_high_edit,'String'));
% Initialize variable determining the number of bins used for histograms
handles.hist_bins = str2double(get(handles.hist_bin_edit,'String'));
% Initialize L_min variable
handles.L_min = str2double(get(handles.L_min_edit,'String'));
% Initialize L_max variable
handles.L_max = str2double(get(handles.L_max_edit,'String'));
% Initialize the frequency resolution variable
handles.frequency_resolution = ...
    str2double(get(handles.frequency_resolution_edit,'String'));
% Initialize inclusion of complex results to off and set respective
% checkbox to the assigend off-value
handles.include_complex_results = 0;
set(handles.include_complex_check,'Value',...
    handles.include_complex_results);
% Switch to not apply the Tracedropper rejection of filaments
handles.non_dropped = false; % If set to true, dropped filaments are included


% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ResultBundleManager wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ResultBundleManager_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in bundle_from_directory_button.
function bundle_from_directory_button_Callback(hObject, eventdata, handles)
% hObject    handle to bundle_from_directory_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% ----
% Get directory that should be turned into a result bundle

% Check if path has already been chosen before, otherwise assign with
% current directory
if handles.source_path == false
    handles.source_path = pwd;
end

% Choose path to directory that should be bundled
handles.source_path = uigetdir(handles.source_path, ...
    'Choose root directory to turn into bundle.');

% ----
% Get directory to save the result bundle to

% Check if target path has already been chosen before, otherwise assign
% with source path
if handles.bundle_path == false
    handles.bundle_path = handles.source_path;
end

[save_file,save_path] = uiputfile(handles.bundle_path, ...
    'Where to save? Save as *.mat file!');
handles.bundle_path = [save_path save_file];
% Make sure that file is saved as a .mat file
if ~strcmp(handles.bundle_path(end-3:end),'.mat')
    handles.bundle_path = [handles.bundle_path '.mat'];
end

% ----
% Create bundle from a directory
directory_to_bundle(handles.source_path,handles.bundle_path);

% ----
% Load the bundle just saved into the manager, the function call makes sure
% that all processes associated with a proper bundle load are executed
handles = load_bundle(hObject,handles);

guidata(hObject,handles)


% --- Executes on button press in select_bundle_button.
function select_bundle_button_Callback(hObject, eventdata, handles)
% hObject    handle to select_bundle_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Check if path has already been chosen before, otherwise assign with
% current directory
if ~handles.bundle_path
    handles.bundle_path = pwd;
end

% Set the new bundle path from user dialog to pick file
[load_file,load_path] = uigetfile([handles.bundle_path], ...
    'Select bundle file to load.');
handles.bundle_path = [load_path load_file];

% Calll to function that executes things associated with loading of a
% bundle
handles = load_bundle(hObject,handles);

guidata(hObject,handles)



function bundle_path_edit_Callback(hObject, eventdata, handles)
% hObject    handle to bundle_path_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bundle_path_edit as text
%        str2double(get(hObject,'String')) returns contents of bundle_path_edit as a double

set(handles.bundle_path_edit,'String',handles.bundle_path)

guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function bundle_path_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bundle_path_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in VL_diagram_button.
function VL_diagram_button_Callback(hObject, eventdata, handles)
% hObject    handle to VL_diagram_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA) 


% Call to function that updates the three keyword lists from the GUI
handles = poll_keyword_edits(hObject,handles);

% Pull trace velocities out of the queried result sections
pull_velocity = ...
    @(result_section) [result_section.trace_results.trace_velocity];
trace_velocity = extract_by_keywords(handles.bundle_path, ...
    handles.any_keywords,handles.with_keywords,handles.not_keywords, ...
    pull_velocity,handles.non_dropped);
trace_velocity = [trace_velocity{:}];

% Pull filament average lenghts out of the queried result sections
pull_average_length = ...
    @(result_section) ...
    [result_section.trace_results.average_filament_length];
average_filament_length = extract_by_keywords(handles.bundle_path, ...
    handles.any_keywords,handles.with_keywords,handles.not_keywords, ...
    pull_average_length,handles.non_dropped);
average_filament_length = [average_filament_length{:}];

% Apply the L_min and L_max limits
in_limits_inds = find( ...
    average_filament_length >= handles.L_min ...
    & average_filament_length <= handles.L_max);

% Find indices for which the filament length has no imaginary part, i.e.
% the indices of all filaments for which the rectangular transformation has
% a real solution
non_imag_inds = in_limits_inds(...
    imag(average_filament_length(in_limits_inds))==0);

% Display average filament length on command line
fprintf('Accurate length average:\n')
fprintf('%f+-%f(StDev), n=%d\n',...
    mean(average_filament_length(non_imag_inds)),...
    std(average_filament_length(non_imag_inds)),...
    numel(non_imag_inds))
fprintf('Qualitative length average:\n')
fprintf('%f+-%f(StDev), n=%d\n\n',...
    mean(average_filament_length(in_limits_inds)),...
    std(average_filament_length(in_limits_inds)),...
    numel(in_limits_inds))


% Plot the (V,L) diagram
VL_figure = figure;
VL_axes = axes('Parent',VL_figure);
plot(VL_axes, ...
    real(average_filament_length), ...
    real(trace_velocity), ...
    'ko','MarkerSize',4)
set(VL_axes,'NextPlot','Add')
plot(VL_axes, ...
    real(average_filament_length(in_limits_inds)), ...
    real(trace_velocity(in_limits_inds)), ...
    'r+','MarkerSize',4,'MarkerFaceColor','none')
plot(VL_axes, ...
    real(average_filament_length(non_imag_inds)), ...
    real(trace_velocity(non_imag_inds)), ...
    'bs','MarkerSize',8,'MarkerFaceColor','none')
plot(VL_axes,handles.L_min.*[1 1], ...
    [0 1.05.*max(trace_velocity)], ...
    'r-')
if isfinite(handles.L_max)
    plot(VL_axes,handles.L_max.*[1 1], ...
        [0 1.05.*max(trace_velocity)], ...
        'r-')
end
set(VL_axes,'YLim', ...
    [0 1.05.*max(trace_velocity)])
set(VL_axes,'NextPlot','Replace')
title(VL_axes,sprintf('v_{max} = %f +- %f (SDev), n_{trcs}=%d', ...
    mean(trace_velocity(in_limits_inds)), ...
    std(trace_velocity(in_limits_inds)), ...
    numel(trace_velocity(in_limits_inds))),...
    'FontSize',12)
xlabel(VL_axes,'Filament length L[\mum]','FontSize',12)
ylabel(VL_axes,'Trace velocity V[\mum/s]','FontSize',12)


function handles = load_bundle(hObject,handles)

fprintf('Started loading result bundle into manager...')

% Set the new bundle path in the GUI
set(handles.bundle_path_edit,'String',handles.bundle_path)

% Get keyword table from selected bundle
keyword_table_load = load(handles.bundle_path,'keyword_table');
handles.keyword_table = keyword_table_load.keyword_table;

% Pull filament average lenghts out of the queried result sections

% Number of sections in the bundle
handles.sections_in_bundle = numel(handles.keyword_table)./2;
% Check if all sections in the bundle have been treated by the TraceDropper
tracedrop_mask_isfield_function = ...
    @(section) isfield(section,{'tracedrop_mask'});
non_dropped_true = true; % To prevent conflict when extracting results here
handles.tracedropped_sections_inds = ...
    extract_by_keywords(handles.bundle_path, ...
    handles.any_keywords,handles.with_keywords,handles.not_keywords, ...
    tracedrop_mask_isfield_function,non_dropped_true);
handles.tracedropped_sections_inds = ...
    [handles.tracedropped_sections_inds{:}];

if sum(handles.tracedropped_sections_inds) == ...
        handles.sections_in_bundle
    set(handles.bundle_status_text,'String',...
        sprintf('%d sections loaded.',handles.sections_in_bundle),...
        'BackgroundColor', ...
        get(0,'defaultUicontrolBackgroundColor'))
else
    set(handles.bundle_status_text,'String',...
        sprintf(['%d sections loaded,' ...
        'not all sections trace-dropped.'],handles.sections_in_bundle),...
        'BackgroundColor', ...
        [1,0,0])
end

fprintf('done.\n')

function handles = poll_keyword_edits(hObject,handles)


% Check if any_keywords have been edited
if handles.any_keywords_edited
    % Get the any_keywords
    any_keywords_string = get(handles.any_keywords_edit,'String');
    if isempty(any_keywords_string)
        handles.any_keywords = {};
    else
        any_keywords_number = sum(any_keywords_string==',')+1;
        any_keywords = textscan( any_keywords_string, '%s' , ...
            any_keywords_number, ...
            'Delimiter',',');
        handles.any_keywords = any_keywords{:};
    end
end

% Check if with_keywords have been edited
if handles.with_keywords_edited
    % Get the with_keywords
    with_keywords_string = get(handles.with_keywords_edit,'String');
    if isempty(with_keywords_string)
        handles.with_keywords = {};
    else
        with_keywords_number = sum(with_keywords_string==',')+1;
        with_keywords = textscan( with_keywords_string, '%s' , ...
            with_keywords_number, ...
            'Delimiter',',');
        handles.with_keywords = with_keywords{:};
    end
end

% Check if not_keywords have been edited
if handles.not_keywords_edited
    % Get the not_keywords
    not_keywords_string = get(handles.not_keywords_edit,'String');
    if isempty(not_keywords_string)
        handles.not_keywords = {};
    else
        not_keywords_number = sum(not_keywords_string==',')+1;
        not_keywords = textscan( not_keywords_string, '%s' , ...
            not_keywords_number, ...
            'Delimiter',',');
        handles.not_keywords = not_keywords{:};
    end
end



function with_keywords_edit_Callback(hObject, eventdata, handles)
% hObject    handle to with_keywords_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.with_keywords_edited = true;

% Hints: get(hObject,'String') returns contents of with_keywords_edit as text
%        str2double(get(hObject,'String')) returns contents of with_keywords_edit as a double

guidata(hObject,handles)



% --- Executes during object creation, after setting all properties.
function with_keywords_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to with_keywords_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function not_keywords_edit_Callback(hObject, eventdata, handles)
% hObject    handle to not_keywords_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.not_keywords_edited = true;

% Hints: get(hObject,'String') returns contents of not_keywords_edit as text
%        str2double(get(hObject,'String')) returns contents of not_keywords_edit as a double

guidata(hObject,handles)



% --- Executes during object creation, after setting all properties.
function not_keywords_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to not_keywords_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function handles = get_extracted_bundle(hObject,handles)
% get_extracted_bundle gets a bundle of result sections that are conform
% with the entered keywords, and contains only sections that have been
% trace-dropped. The extracted bundle is stored in handles.extracted_bundle

% Call to the poll_keyword_edits function to get keywords from GUI edit
% fields
handles = poll_keyword_edits(hObject,handles);

% Call to function that returns a bundle only containing sections that
% confirm with keywords and that have been trace-dropped
handles.extracted_bundle = ...
    extract_by_keywords( handles.bundle, ...
    handles.any_keywords, ...
    handles.with_keywords, ...
    handles.not_keywords );


% --- Executes on button press in Save_extracted_bundle.
function Save_extracted_bundle_Callback(hObject, eventdata, handles)
% hObject    handle to Save_extracted_bundle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isfield(handles,'bundle_path')
    %If no result bundle file has been loaded yet
    msgbox('You must choose or construct a result bundle, see above buttons.', ...
        'No result bundle chosen yet.',...
        'Icon','warn')
else
    % Get the extracted bundle according to keywords, and create a save_bundle
    % from it, which can then be saved
    
    % Update the keywords from keyword entry fields
    handles = poll_keyword_edits(hObject,handles);
    
    % Let the user pick a target file for the extracted bundle
    if isfield(handles,'extracted_save_path')
        [save_file,save_path] = uiputfile(handles.extracted_save_path, ...
            'Where to save the extracted bundle file?');
    else
        [save_file,save_path] = uiputfile(handles.bundle_path, ...
            'Where to save the extracted bundle file?');
    end
    handles.extracted_save_path = [save_path save_file];
    
    % Extract and save a bundle according to current keywords
    bundle_to_queried_bundle( ...
        handles.bundle_path,handles.extracted_save_path, ...
        handles.any_keywords,handles.with_keywords,handles.not_keywords);
    
end

guidata(hObject,handles)


% --- Executes on button press in F2F_hist_button.
function F2F_hist_button_Callback(hObject, eventdata, handles)
% hObject    handle to F2F_hist_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% ----
% Draw frame-to-frame (F2F) velocity histogram for current extracted bundle

% Call to function that updates the three keyword lists from the GUI
handles = poll_keyword_edits(hObject,handles);

% Pull trace velocities out of the queried result sections
pull_velocity = ...
    @(result_section) [result_section.trace_results.trace_velocity];
trace_velocity = extract_by_keywords(handles.bundle_path, ...
    handles.any_keywords,handles.with_keywords,handles.not_keywords, ...
    pull_velocity,handles.non_dropped);
trace_velocity = [trace_velocity{:}];

% Pull filament average lenghts out of the queried result sections
pull_average_length = ...
    @(result_section) ...
    [result_section.trace_results.average_filament_length];
average_filament_length = extract_by_keywords(handles.bundle_path, ...
    handles.any_keywords,handles.with_keywords,handles.not_keywords, ...
    pull_average_length,handles.non_dropped);
average_filament_length = [average_filament_length{:}];

% Apply the L_min and L_max limits
in_limits_inds = find( ...
    average_filament_length >= handles.L_min ...
    & average_filament_length <= handles.L_max);

% Find indices for which the filament length has no imaginary part, i.e.
% the indices of all filaments for which the rectangular transformation has
% a real solution
non_imag_inds = in_limits_inds(...
    imag(average_filament_length(in_limits_inds))==0);

% Pull frame-to-frame velocities out of the queried result sections
pull_f2f_velocities = ...
    @(result_section) ...
    {result_section.trace_results.frame_to_frame_velocities};
per_frame_velocities = extract_by_keywords(handles.bundle_path, ...
    handles.any_keywords,handles.with_keywords,handles.not_keywords, ...
    pull_f2f_velocities,handles.non_dropped);
per_frame_velocities = ...
    [per_frame_velocities{:}];

if handles.include_complex_results
    per_frame_velocities = ...
        [per_frame_velocities{in_limits_inds}];
else
    per_frame_velocities = ...
        [per_frame_velocities{non_imag_inds}];
end

% Make a histogram of the frame-to-frame velocities
binning_edges = linspace(0,handles.V_high,handles.hist_bins+1);
bin_centers = (binning_edges(1:end-1)+binning_edges(2:end))./2;
bin_counts = histc(per_frame_velocities,binning_edges);
hist_figure = figure;
hist_axes = axes('Parent',hist_figure);
hist_handle = bar(bin_centers,bin_counts(1:end-1), .8, ...
    'Parent',hist_axes,...
    'LineStyle','none');
xlabel('F2F V[\mum/s]')
ylabel('Bin count')
set(hist_axes,'NextPlot','Add')

% Add two Gaussian mixture model to histograms
[ proportions,means,sigmas,mixture_model ] = ...
    two_gaussian_fit( ...
    per_frame_velocities,handles.V_high );
vv_support = linspace(0,handles.V_high,500).';
plot(hist_axes,vv_support, ...
    mixture_model(vv_support)./max(mixture_model(vv_support)) ...
    .*max(bin_counts(1:end-1)),...
    'r-')

% Add v_max and motile fraction to high velocitiy peak
text( 1.1.*means(2), ...
    mixture_model(means(2))./max(mixture_model(vv_support)) ...
    .*max(bin_counts(1:end-1)),...
    sprintf('v_{max}=%f\nf_{mot}=%f',means(2),proportions(2)), ...
    'Parent',hist_axes)

title(hist_axes, sprintf('F2F velocities, L in [%f,%f], n=%d', ...
    handles.L_min,handles.L_max, ...
    numel(per_frame_velocities)))
set(hist_axes,'YLim',[0 1.1.*max(bin_counts(1:end-1))])
set(hist_axes,'XLim',[0 handles.V_high])
set(hist_axes,'NextPlot','Replace')



function hist_bin_edit_Callback(hObject, eventdata, handles)
% hObject    handle to hist_bin_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of hist_bin_edit as text
%        str2double(get(hObject,'String')) returns contents of hist_bin_edit as a double

% Make sure that the histogram count is always integer and greater than 1
handles.hist_bins = str2double(get(handles.hist_bin_edit,'String'));
handles.hist_bins = ceil(handles.hist_bins).*(handles.hist_bins>=2) ...
    + 2.*(handles.hist_bins<2);
set(handles.hist_bin_edit,'String',num2str(handles.hist_bins))

guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function hist_bin_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to hist_bin_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function V_high_edit_Callback(hObject, eventdata, handles)
% hObject    handle to V_high_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of V_high_edit as text
%        str2double(get(hObject,'String')) returns contents of V_high_edit as a double

% Make sure that V_high is always positive
% Make sure that the histogram count is always integer and greater than 1
handles.old_V_high = handles.V_high;

handles.V_high = str2double(get(handles.V_high_edit,'String'));
handles.V_high = handles.V_high.*(handles.V_high>0) + ...
    handles.old_V_high.*(handles.V_high<=0);
set(handles.V_high_edit,'String',num2str(handles.V_high))

guidata(hObject,handles)




% --- Executes during object creation, after setting all properties.
function V_high_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to V_high_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Trace_hist_button.
function Trace_hist_button_Callback(hObject, eventdata, handles)
% hObject    handle to Trace_hist_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% ----
% Draw frame-to-frame (F2F) velocity histogram for current extracted bundle

% Call to function that updates the three keyword lists from the GUI
handles = poll_keyword_edits(hObject,handles);

% Pull trace velocities out of the queried result sections
pull_velocity = ...
    @(result_section) [result_section.trace_results.trace_velocity];
trace_velocity = extract_by_keywords(handles.bundle_path, ...
    handles.any_keywords,handles.with_keywords,handles.not_keywords, ...
    pull_velocity,handles.non_dropped);
trace_velocity = [trace_velocity{:}];

% Pull filament average lenghts out of the queried result sections
pull_average_length = ...
    @(result_section) ...
    [result_section.trace_results.average_filament_length];
average_filament_length = extract_by_keywords(handles.bundle_path, ...
    handles.any_keywords,handles.with_keywords,handles.not_keywords, ...
    pull_average_length,handles.non_dropped);
average_filament_length = [average_filament_length{:}];

% Apply the L_min and L_max limits
in_limits_inds = find( ...
    average_filament_length >= handles.L_min ...
    & average_filament_length <= handles.L_max);

% Find indices for which the filament length has no imaginary part, i.e.
% the indices of all filaments for which the rectangular transformation has
% a real solution
non_imag_inds = in_limits_inds(...
    imag(average_filament_length(in_limits_inds))==0);

if handles.include_complex_results
    trace_velocity = ...
        [trace_velocity(in_limits_inds)];
else
    trace_velocities = ...
        [trace_velocity(non_imag_inds)];
end

% Make a histogram of the frame-to-frame velocities
binning_edges = linspace(0,handles.V_high,handles.hist_bins+1);
bin_centers = (binning_edges(1:end-1)+binning_edges(2:end))./2;
bin_counts = histc(trace_velocities,binning_edges);
hist_figure = figure;
hist_axes = axes('Parent',hist_figure);
hist_handle = bar(bin_centers,bin_counts(1:end-1), .8, ...
    'Parent',hist_axes,...
    'LineStyle','none');
xlabel(hist_axes,'Trace V[\mum/s]')
ylabel(hist_axes,'Bin count')
set(hist_axes,'NextPlot','Add')

% Add two Gaussian mixture model to histograms
[ proportions,means,sigmas,mixture_model ] = ...
    two_gaussian_fit( trace_velocities,handles.V_high );
vv_support = linspace(0,handles.V_high,500).';
plot(hist_axes,vv_support, ...
    mixture_model(vv_support)./max(mixture_model(vv_support)) ...
    .*max(bin_counts(1:end-1)),...
    'r-')

% Add v_max and motile fraction to high velocitiy peak
text( 1.1.*means(2), ...
    mixture_model(means(2))./max(mixture_model(vv_support)) ...
    .*max(bin_counts(1:end-1)),...
    sprintf('v_{max}=%f\nf_{mot}=%f',means(2),proportions(2)), ...
    'Parent',hist_axes)

title(hist_axes, sprintf('Trace velocities, L in [%f,%f], n=%d', ...
    handles.L_min,handles.L_max, ...
    numel(trace_velocity)))
set(hist_axes,'YLim',[0 1.1.*max(bin_counts(1:end-1))])
set(hist_axes,'XLim',[0 handles.V_high])
set(hist_axes,'NextPlot','Replace')

guidata(hObject,handles)





function L_min_edit_Callback(hObject, eventdata, handles)
% hObject    handle to L_min_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of L_min_edit as text
%        str2double(get(hObject,'String')) returns contents of L_min_edit as a double

% Make sure that L_min is always positive
handles.L_min = str2double(get(handles.L_min_edit,'String'));
handles.L_min = handles.L_min.*(handles.L_min>0);
set(handles.L_min_edit,'String',num2str(handles.L_min))

guidata(hObject,handles)



% --- Executes during object creation, after setting all properties.
function L_min_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to L_min_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function L_max_edit_Callback(hObject, eventdata, handles)
% hObject    handle to L_max_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of L_max_edit as text
%        str2double(get(hObject,'String')) returns contents of L_max_edit as a double

% Make sure that L_max is always greater than L_min
handles.L_max = str2double(get(handles.L_max_edit,'String'));
handles.L_max = handles.L_max.*(handles.L_max>handles.L_min);
if handles.L_max <= handles.L_min
    handles.L_max = Inf;
end
set(handles.L_max_edit,'String',num2str(handles.L_max))

guidata(hObject,handles)



% --- Executes during object creation, after setting all properties.
function L_max_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to L_max_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in frequency_button.
function frequency_button_Callback(hObject, eventdata, handles)
% hObject    handle to frequency_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% ----
% Draw frequency spectrum for current extracted bundle

% Call to function that updates the three keyword lists from the GUI
handles = poll_keyword_edits(hObject,handles);

% Pull trace velocities out of the queried result sections
pull_velocity = ...
    @(result_section) [result_section.trace_results.trace_velocity];
trace_velocity = extract_by_keywords(handles.bundle_path, ...
    handles.any_keywords,handles.with_keywords,handles.not_keywords, ...
    pull_velocity,handles.non_dropped);
trace_velocity = [trace_velocity{:}];

% Pull filament average lenghts out of the queried result sections
pull_average_length = ...
    @(result_section) ...
    [result_section.trace_results.average_filament_length];
average_filament_length = extract_by_keywords(handles.bundle_path, ...
    handles.any_keywords,handles.with_keywords,handles.not_keywords, ...
    pull_average_length,handles.non_dropped);
average_filament_length = [average_filament_length{:}];

% Apply the L_min and L_max limits
in_limits_inds = find( ...
    average_filament_length >= handles.L_min ...
    & average_filament_length <= handles.L_max);

% Find indices for which the filament length has no imaginary part, i.e.
% the indices of all filaments for which the rectangular transformation has
% a real solution
non_imag_inds = in_limits_inds(...
    imag(average_filament_length(in_limits_inds))==0);

% Pull frame-to-frame velocities out of the queried result sections
pull_f2f_velocities = ...
    @(result_section) ...
    {result_section.trace_results.frame_to_frame_velocities};
per_frame_velocities = extract_by_keywords(handles.bundle_path, ...
    handles.any_keywords,handles.with_keywords,handles.not_keywords, ...
    pull_f2f_velocities,handles.non_dropped);
per_frame_velocities = [per_frame_velocities{:}];

if handles.include_complex_results
    per_frame_velocities = ...
        per_frame_velocities(in_limits_inds);
else
    per_frame_velocities = ...
        per_frame_velocities(non_imag_inds);
end

% Draw figure with spectra of all traces
pull_frame_rate = @(section) section.video_properties.merged_frame_rate;
fprintf(['Frame rates of different sections,\n' ...
    'Bad if they are different!\n'])
frame_rate = extract_by_keywords(handles.bundle_path, ...
    handles.any_keywords,handles.with_keywords,handles.not_keywords, ...
    pull_frame_rate,handles.non_dropped);
frame_rate = [frame_rate{:}];
disp(frame_rate)
frame_rate = mean(frame_rate(:));

% Extract frequency spectrum using Welch's method
window_width = 2.*(handles.frequency_resolution-1);
window_overlap = []; fft_points = window_width;
sampling_frequency = frame_rate;
power_spectrum_array = [];

long_enough_velocity_traces = 0;
for kk = 1:numel(per_frame_velocities)
    this_velocities = ...
        real(per_frame_velocities{kk});
    if numel(this_velocities)>window_width
        long_enough_velocity_traces = long_enough_velocity_traces+1;
        [power_spectrum,frequencies] = ...
            pwelch(this_velocities,window_width, ...
            window_overlap,fft_points,sampling_frequency);
        power_spectrum_array = [power_spectrum_array;power_spectrum.'];
    end
end

if long_enough_velocity_traces > 0
    % If at least one filament could be analyzed for its frequency spectrum
    
    power_spectrum_mean = mean(power_spectrum_array,1);
    power_spectrum_stdev = std(power_spectrum_array,1);
    
    power_to_dezibel = @(power_value) 10.*log10(power_value);
    
    spectrum_figure = figure;
    spectrum_axes = axes('Parent',spectrum_figure);
    
    set(spectrum_axes,'NextPlot','Add')
    plot(frequencies,power_to_dezibel(power_spectrum_mean),'ko-',...
        'LineWidth',1.5)
    plot(frequencies, ...
        power_to_dezibel(power_spectrum_mean+power_spectrum_stdev),'r-',...
        'LineWidth',1)
    xlabel(spectrum_axes,'Frequency[1/s]','FontSize',12)
    ylabel(spectrum_axes,'Power[dB]','FontSize',12)
    title(spectrum_axes,sprintf('Considered traces: %d',...
        long_enough_velocity_traces))
    set(spectrum_axes,'NextPlot','Replace', ...
        'Box','on',...
        'YScale','linear')
else
    msgbox(['No traces fulfill criteria for frequency analysis.' ...
        'Try reducing frequency resolution.'])
end

guidata(hObject,handles)



function frequency_resolution_edit_Callback(hObject, eventdata, handles)
% hObject    handle to frequency_resolution_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of frequency_resolution_edit as text
%        str2double(get(hObject,'String')) returns contents of frequency_resolution_edit as a double

% Make sure that the frequency resolution is always 2 or higher

handles.frequency_resolution = ...
    str2double(get(handles.frequency_resolution_edit,'String'));
handles.frequency_resolution = ...
    (handles.frequency_resolution>1).*ceil(handles.frequency_resolution) ...
    + 2.*(handles.frequency_resolution<=1);
set(handles.frequency_resolution_edit,'String', ...
    num2str(handles.frequency_resolution))

guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function frequency_resolution_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to frequency_resolution_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in SL_diagram_button.
function SL_diagram_button_Callback(hObject, eventdata, handles)
% hObject    handle to SL_diagram_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Call to function that updates the three keyword lists from the GUI
handles = poll_keyword_edits(hObject,handles);

% Pull trace velocities out of the queried result sections
pull_velocity = ...
    @(result_section) [result_section.trace_results.trace_velocity];
trace_velocity = extract_by_keywords(handles.bundle_path, ...
    handles.any_keywords,handles.with_keywords,handles.not_keywords, ...
    pull_velocity,handles.non_dropped);
trace_velocity = [trace_velocity{:}];

% Pull filament average lenghts out of the queried result sections
pull_average_length = ...
    @(result_section) ...
    [result_section.trace_results.average_filament_length];
average_filament_length = extract_by_keywords(handles.bundle_path, ...
    handles.any_keywords,handles.with_keywords,handles.not_keywords, ...
    pull_average_length,handles.non_dropped);
average_filament_length = [average_filament_length{:}];

% Apply the L_min and L_max limits
in_limits_inds = find( ...
    average_filament_length >= handles.L_min ...
    & average_filament_length <= handles.L_max);

% Find indices for which the filament length has no imaginary part, i.e.
% the indices of all filaments for which the rectangular transformation has
% a real solution
non_imag_inds = in_limits_inds(...
    imag(average_filament_length(in_limits_inds))==0);

% Pull frame-to-frame velocities out of the queried result sections
pull_f2f_velocities = ...
    @(result_section) ...
    {result_section.trace_results.frame_to_frame_velocities};
per_frame_velocities = extract_by_keywords(handles.bundle_path, ...
    handles.any_keywords,handles.with_keywords,handles.not_keywords, ...
    pull_f2f_velocities,handles.non_dropped);
per_frame_velocities = [per_frame_velocities{:}];

% Pull frame and mean solidities out of the queried result sections
pull_frame_solidities = ...
    @(result_section) ...
    {result_section.trace_results.filament_solidities};
frame_solidities = extract_by_keywords(handles.bundle_path, ...
    handles.any_keywords,handles.with_keywords,handles.not_keywords, ...
    pull_frame_solidities,handles.non_dropped);
frame_solidities = ...
    [frame_solidities{:}];
mean_solidity = ...
    cellfun(@(xx) mean(xx(:)),frame_solidities);

% Pull trace solidities out of the queried result sections
pull_trace_solidity = ...
    @(result_section) ...
    [result_section.trace_results.trace_solidity];
trace_solidity = extract_by_keywords(handles.bundle_path, ...
    handles.any_keywords,handles.with_keywords,handles.not_keywords, ...
    pull_trace_solidity,handles.non_dropped);
trace_solidity = ...
    [trace_solidity{:}];




% if handles.include_complex_results
%     per_frame_velocities = ...
%         per_frame_velocities(in_limits_inds);
%     trace_velocity = trace_velocity(in_limits_inds);
%     average_filament_length = average_filament_length(in_limits_inds);
%     mean_solidity = mean_solidity(in_limits_inds);
%     trace_solidity = trace_solidity(in_limits_inds);
% else
%     per_frame_velocities = ...
%         per_frame_velocities(non_imag_inds);
%     trace_velocity = trace_velocity(non_imag_inds);
%     average_filament_length = average_filament_length(non_imag_inds);
%     mean_solidity = mean_solidity(non_imag_inds);
%     trace_solidity = trace_solidity(non_imag_inds);
% end




% Plot the (S,L) diagram, mean solidity over individual frames
SL_figure = figure;
SL_mean_axes = subplot(3,1,1,'Parent',SL_figure);
plot(SL_mean_axes, ...
    real(average_filament_length), ...
    real(mean_solidity), ...
    'ko','MarkerSize',4)
set(SL_mean_axes,'NextPlot','Add')
plot(SL_mean_axes, ...
    real(average_filament_length(in_limits_inds)), ...
    real(mean_solidity(in_limits_inds)), ...
    'r+','MarkerSize',4,'MarkerFaceColor','none')
plot(SL_mean_axes, ...
    real(average_filament_length(non_imag_inds)), ...
    real(mean_solidity(non_imag_inds)), ...
    'bs','MarkerSize',8,'MarkerFaceColor','none')
plot(SL_mean_axes,handles.L_min.*[1 1], ...
    [0 1.05.*max(mean_solidity)], ...
    'r-')
if isfinite(handles.L_max)
    plot(SL_mean_axes,handles.L_max.*[1 1], ...
        [0 1.05.*max(mean_solidity)], ...
        'r-')
end
set(SL_mean_axes,'YLim', ...
    [0 1.05.*max(mean_solidity)])
set(SL_mean_axes,'NextPlot','Replace')
title(SL_mean_axes,sprintf('S = %f +- %f (SDev), n_{trcs}=%d', ...
    mean(mean_solidity(non_imag_inds)), ...
    std(mean_solidity(non_imag_inds)), ...
    numel(mean_solidity(non_imag_inds))),...
    'FontSize',12)
xlabel(SL_mean_axes,'Filament length L[\mum]','FontSize',12)
ylabel(SL_mean_axes,'Filament average solidity','FontSize',12)

% Plot the histogram of all individual frame solidities
SL_hist_axes = subplot(3,1,2,'Parent',SL_figure);

% Make a histogram of the frame-to-frame velocities
binning_edges = linspace(0,1,handles.hist_bins+1);
bin_centers = (binning_edges(1:end-1)+binning_edges(2:end))./2;
pooled_solidities = [frame_solidities{non_imag_inds}];
bin_counts = histc(pooled_solidities,binning_edges);
max_bin_count = max(bin_counts);
hist_handle = bar(bin_centers,bin_counts(1:end-1), .8, ...
    'Parent',SL_hist_axes,...
    'LineStyle','none');
xlabel(SL_hist_axes,'Frame solidity','FontSize',12)
ylabel(SL_hist_axes,'Bin count','FontSize',12)
mean_solidity = mean(pooled_solidities);
stddev_solidity = std(pooled_solidities);
title(SL_hist_axes,sprintf('Mean solidity: %f+-%f(SDev,n=%d)',...
    mean_solidity,stddev_solidity,numel(pooled_solidities)),...
    'FontSize',12)
set(SL_hist_axes,'NextPlot','Add',...
    'YLim',[0 max_bin_count.*1.1])
plot(mean_solidity.*[1 1],[0 max_bin_count.*1.1],'k--','LineWidth',1.25)


% Plot the (S,L) diagram, solidity of whole traces
SL_trace_axes = subplot(3,1,3,'Parent',SL_figure);
plot(SL_trace_axes, ...
    real(average_filament_length), ...
    real(trace_solidity), ...
    'ko','MarkerSize',4)
set(SL_trace_axes,'NextPlot','Add')
plot(SL_trace_axes, ...
    real(average_filament_length(in_limits_inds)), ...
    real(trace_solidity(in_limits_inds)), ...
    'r+','MarkerSize',4,'MarkerFaceColor','none')
plot(SL_trace_axes, ...
    real(average_filament_length(non_imag_inds)), ...
    real(trace_solidity(non_imag_inds)), ...
    'bs','MarkerSize',8,'MarkerFaceColor','none')
plot(SL_trace_axes,handles.L_min.*[1 1], ...
    [0 1.05.*max(trace_solidity)], ...
    'r-')
if isfinite(handles.L_max)
    plot(SL_trace_axes,handles.L_max.*[1 1], ...
        [0 1.05.*max(trace_solidity)], ...
        'r-')
end
set(SL_trace_axes,'YLim', ...
    [0 1.05.*max(trace_solidity)])
set(SL_trace_axes,'NextPlot','Replace')
title(SL_trace_axes,sprintf('S = %f +- %f (SDev), n_{trcs}=%d', ...
    mean(trace_solidity(non_imag_inds)), ...
    std(trace_solidity(non_imag_inds)), ...
    numel(trace_solidity(non_imag_inds))),...
    'FontSize',12)
xlabel(SL_trace_axes,'Filament length L[\mum]','FontSize',12)
ylabel(SL_trace_axes,'Trace solidity','FontSize',12)

guidata(hObject,handles)


% --- Executes on button press in pixel_vs_centroid_button.
function pixel_vs_centroid_button_Callback(hObject, eventdata, handles)
% hObject    handle to pixel_vs_centroid_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Collect data for and then construct and plot a pixel displacement vs.
% centroid displacement plot

% Call to function that updates the three keyword lists from the GUI
handles = poll_keyword_edits(hObject,handles);

% Pull trace velocities out of the queried result sections
pull_velocity = ...
    @(result_section) [result_section.trace_results.trace_velocity];
trace_velocity = extract_by_keywords(handles.bundle_path, ...
    handles.any_keywords,handles.with_keywords,handles.not_keywords, ...
    pull_velocity,handles.non_dropped);
trace_velocity = [trace_velocity{:}];

% Pull filament average lenghts out of the queried result sections
pull_average_length = ...
    @(result_section) ...
    [result_section.trace_results.average_filament_length];
average_filament_length = extract_by_keywords(handles.bundle_path, ...
    handles.any_keywords,handles.with_keywords,handles.not_keywords, ...
    pull_average_length,handles.non_dropped);
average_filament_length = [average_filament_length{:}];

% Apply the L_min and L_max limits
in_limits_inds = find( ...
    average_filament_length >= handles.L_min ...
    & average_filament_length <= handles.L_max);

% Find indices for which the filament length has no imaginary part, i.e.
% the indices of all filaments for which the rectangular transformation has
% a real solution
non_imag_inds = in_limits_inds(...
    imag(average_filament_length(in_limits_inds))==0);

% Pull frame-to-frame velocities out of the queried result sections
pull_f2f_velocities = ...
    @(result_section) ...
    {result_section.trace_results.frame_to_frame_velocities};
per_frame_velocities = extract_by_keywords(handles.bundle_path, ...
    handles.any_keywords,handles.with_keywords,handles.not_keywords, ...
    pull_f2f_velocities,handles.non_dropped);
per_frame_velocities = [per_frame_velocities{:}];

% Pull frame and mean solidities out of the queried result sections
pull_frame_solidities = ...
    @(result_section) ...
    {result_section.trace_results.filament_solidities};
frame_solidities = extract_by_keywords(handles.bundle_path, ...
    handles.any_keywords,handles.with_keywords,handles.not_keywords, ...
    pull_frame_solidities,handles.non_dropped);
frame_solidities = ...
    [frame_solidities{:}];
mean_solidity = ...
    cellfun(@(xx) mean(xx(:)),frame_solidities);

% Pull trace solidities out of the queried result sections
pull_trace_solidity = ...
    @(result_section) ...
    [result_section.trace_results.trace_solidity];
trace_solidity = extract_by_keywords(handles.bundle_path, ...
    handles.any_keywords,handles.with_keywords,handles.not_keywords, ...
    pull_trace_solidity,handles.non_dropped);
trace_solidity = ...
    [trace_solidity{:}];

% Function to extract pixel displacements from cell array containing the
% filament pixels in each of the frames
get_pixel_displacements = @(frame_ind_cell) ...
    cellfun(@(old_inds,new_inds) numel(setdiff(new_inds,old_inds)), ...
    frame_ind_cell(1:end-1),frame_ind_cell(2:end));
% Pull pixel displacements using above function
pull_displacements = ...
    @(result_section) cellfun(get_pixel_displacements,...
    {result_section.trace_results.filament_pixel_indices_linear},...
    'UniformOutput',false);
pixel_displacements = ...
    extract_by_keywords(handles.bundle_path, ...
    handles.any_keywords,handles.with_keywords,handles.not_keywords,...
    pull_displacements,handles.non_dropped);
pixel_displacements = ...
    [pixel_displacements{:}];

% Sort by ascending length
[average_filament_length,sort_inds] = ...
    sort(real(average_filament_length),'ascend');
per_frame_velocities = ...
    per_frame_velocities(sort_inds);
pixel_displacements = ...
    pixel_displacements(sort_inds);

% Construct pixel displacement over centroid frame-to-frame ratio for each
% individual frame for each individual filament
ctrd_over_pxl = cellfun(@(ctrd,pxl) ctrd./pxl,...
    per_frame_velocities,...
    pixel_displacements,...
    'UniformOutput',false); % Get the ratio
for kk = 1:numel(ctrd_over_pxl)
    %Set non-defined 0/0 to zero
    ctrd_over_pxl{kk}(~isfinite(ctrd_over_pxl{kk})) = 0;
end
% Get the necessary data to plot average with standard deviation
mean_ctrd_over_pxl = cellfun(@(xx)mean(xx(:)),ctrd_over_pxl);
stddev_ctrd_over_pxl = cellfun(@(xx)std(xx(:)),ctrd_over_pxl);
length_support = cellfun(@(ctrd_over_pxl,avg_lengths)...
    avg_lengths.*ones(size(ctrd_over_pxl)),...
    ctrd_over_pxl,num2cell(average_filament_length),...
    'UniformOutput',false);


% Apply the L_min and L_max limits
in_limits_inds = find( ...
    average_filament_length >= handles.L_min ...
    & average_filament_length <= handles.L_max);

% Find indices for which the filament length has no imaginary part, i.e.
% the indices of all filaments for which the rectangular transformation has
% a real solution
non_imag_inds = in_limits_inds(...
    imag(average_filament_length(in_limits_inds))==0);


% Plot pixel change vs. frame-to-frame velocities
ctrd_vs_pxl_figure = figure;
ctrd_vs_pxl_axes = axes('Parent',ctrd_vs_pxl_figure);
% All data points
plot(ctrd_vs_pxl_axes,...
    [length_support{:}],[ctrd_over_pxl{:}], ...
    'ko','MarkerSize',4,'MarkerEdgeColor',[0.25,0.25,0.25])
set(ctrd_vs_pxl_axes,'NextPlot','Add')
% Only valid data points
plot(ctrd_vs_pxl_axes, ...
    [length_support{non_imag_inds}],[ctrd_over_pxl{non_imag_inds}], ...
    'r+','MarkerSize',4,'MarkerFaceColor','none')
% Mean curve
plot(ctrd_vs_pxl_axes, ...
    real(average_filament_length), ...
    mean_ctrd_over_pxl, ...
    'k-','LineWidth',2)
plot(ctrd_vs_pxl_axes, ...
    real(average_filament_length), ...
    mean_ctrd_over_pxl+stddev_ctrd_over_pxl, ...
    'k--','LineWidth',1)
plot(ctrd_vs_pxl_axes, ...
    real(average_filament_length), ...
    mean_ctrd_over_pxl-stddev_ctrd_over_pxl, ...
    'k--','LineWidth',1)
% Length limit lines
plot(ctrd_vs_pxl_axes,handles.L_min.*[1 1], ...
    [0 1.05.*max([ctrd_over_pxl{:}])], ...
    'r-')
if isfinite(handles.L_max)
    plot(ctrd_vs_pxl_axes,handles.L_max.*[1 1], ...
        [0 1.05.*max([ctrd_over_pxl{:}])], ...
        'r-')
end

set(ctrd_vs_pxl_axes,'YLim', ...
    [0 1.05.*max(mean_ctrd_over_pxl+stddev_ctrd_over_pxl)])
set(ctrd_vs_pxl_axes,'NextPlot','Replace')
title(ctrd_vs_pxl_axes,sprintf('F2F/pxls = %f +- %f (SDev), n_{frms}=%d', ...
    mean([ctrd_over_pxl{non_imag_inds}]), ...
    std([ctrd_over_pxl{non_imag_inds}]), ...
    numel([ctrd_over_pxl{non_imag_inds}])),...
    'FontSize',12)
ylabel(ctrd_vs_pxl_axes,'F2F/pxls[\mum/s]','FontSize',12)
xlabel(ctrd_vs_pxl_axes,'Filament length L[\mum]','FontSize',12)

guidata(hObject,handles)



function any_keywords_edit_Callback(hObject, eventdata, handles)
% hObject    handle to any_keywords_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.any_keywords_edited = true;

% Hints: get(hObject,'String') returns contents of any_keywords_edit as text
%        str2double(get(hObject,'String')) returns contents of any_keywords_edit as a double

guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function any_keywords_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to any_keywords_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in include_complex_check.
function include_complex_check_Callback(hObject, eventdata, handles)
% hObject    handle to include_complex_check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of include_complex_check

handles.include_complex_results = ...
    get(handles.include_complex_check,'Value');

guidata(hObject,handles)


% --- Executes on button press in include_nondropped_check.
function include_nondropped_check_Callback(hObject, eventdata, handles)
% hObject    handle to include_nondropped_check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of include_nondropped_check

handles.non_dropped = ...
    get(handles.include_nondropped_check,'Value');

guidata(hObject,handles)