function varargout = CollectResultFiles(varargin)
% COLLECTRESULTFILES MATLAB code for CollectResultFiles.fig
%      COLLECTRESULTFILES, by itself, creates a new COLLECTRESULTFILES or raises the existing
%      singleton*.
%
%      H = COLLECTRESULTFILES returns the handle to a new COLLECTRESULTFILES or the handle to
%      the existing singleton*.
%
%      COLLECTRESULTFILES('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in COLLECTRESULTFILES.M with the given input arguments.
%
%      COLLECTRESULTFILES('Property','Value',...) creates a new COLLECTRESULTFILES or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before CollectResultFiles_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to CollectResultFiles_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help CollectResultFiles

% Last Modified by GUIDE v2.5 31-Jul-2012 16:03:50

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @CollectResultFiles_OpeningFcn, ...
                   'gui_OutputFcn',  @CollectResultFiles_OutputFcn, ...
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


% --- Executes just before CollectResultFiles is made visible.
function CollectResultFiles_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to CollectResultFiles (see VARARGIN)

% Choose default command line output for CollectResultFiles
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes CollectResultFiles wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = CollectResultFiles_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


function [v_mean,v_STD,L_mean,L_STD,f_mot,nn] = ...
    get_motility_quantifiers(L_min,L_max,V_min,V_max,V_thresh,...
    result_paths,complex_len_in,complex_vel_in,non_dropped,do_plots)

number_of_results = numel(result_paths);

lengths = [];
velocities = [];
f2f_velocities = {};

handles.prog_bar = waitbar(0,'Starting data read in...');

if non_dropped
    % Include all filaments, irrespective of TraceDropper application
    for kk = 1:number_of_results
        
        handles.prog_bar = ...
            waitbar(kk./number_of_results,handles.prog_bar,'Reading in data');
        
        % Extract the trace results for this result file
        this_trace_results = load(result_paths{kk},'trace_results');
        this_trace_results = this_trace_results.trace_results.trace_results;

        this_trace_results = this_trace_results;
        lengths = [lengths,this_trace_results.average_filament_length];
        velocities = [velocities,this_trace_results.trace_velocity];
        f2f_velocities = ...
            [f2f_velocities,{this_trace_results.frame_to_frame_velocities}];
    end
else
    for kk = 1:number_of_results
        % Include only filaments that were found OK in the TraceDropper
        
        handles.prog_bar = ...
            waitbar(kk./number_of_results,handles.prog_bar,'Reading in data');
        
        % Extract the tracedrop-mask for this result file
        tracedrop_mask = load(result_paths{kk},'tracedrop_mask');
        % Check if a tracedrop mask actually existed
        if isfield(tracedrop_mask,'tracedrop_mask')
            tracedrop_mask = tracedrop_mask.tracedrop_mask;
            
            % Extract the trace results for this result file
            this_trace_results = load(result_paths{kk},'trace_results');
            this_trace_results = this_trace_results.trace_results.trace_results;
            
            this_trace_results = this_trace_results(tracedrop_mask);
            lengths = [lengths,this_trace_results.average_filament_length];
            velocities = [velocities,this_trace_results.trace_velocity];
            f2f_velocities = ...
                [f2f_velocities,{this_trace_results.frame_to_frame_velocities}];
            
        else
            fprintf(['Attention: traces were not dropped for file %d,\n',...
                '%s\nAll results from this file are disregarded.\n'],...
                kk,result_paths{kk})
        end
        
    end
end

if do_plots
    figure
    subplot(2,1,1)
    plot(real(lengths),real(velocities),'ko')
    max_length = max(real(lengths));max_velocity = max(real(velocities));
    hold on
end

if complex_len_in && complex_vel_in
    in_inds = real(lengths)>=L_min & real(lengths)<=L_max ...
        & real(velocities)>=V_min & real(velocities)<=V_max;
elseif complex_len_in
    in_inds = ~imag(lengths) ...
        & real(lengths)>=L_min & real(lengths)<=L_max ...
        & real(velocities)>=V_min & real(velocities)<=V_max;
elseif complex_vel_in
    in_inds = ~imag(velocities) ...
        & real(lengths)>=L_min & real(lengths)<=L_max ...
        & real(velocities)>=V_min & real(velocities)<=V_max;
else
    in_inds = ~imag(lengths) & ~imag(velocities) ...
        & real(lengths)>=L_min & real(lengths)<=L_max ...
        & real(velocities)>=V_min & real(velocities)<=V_max;
end

% Pool the basic properties together for mean and standard deviation
% calculation
lengths = lengths(in_inds);
velocities = velocities(in_inds);
% number of filaments
nn = numel(lengths);

if nn>0
    % If at least one filament is left over to calculate results for
    
    % Trace velocity quantifiers
    v_mean = mean(real(velocities)); v_STD = std(real(velocities));
    L_mean = mean(real(lengths)); L_STD = std(real(lengths));
    
    % Motile fraction
    f2f_velocities = [f2f_velocities{in_inds}];
    f_mot = sum(f2f_velocities>=V_thresh)./numel(f2f_velocities);
    
    if do_plots
        xlabel('L[\mum]')
        ylabel('V[\mum/s]')
        set(gca,'XLim',[0 max_length],'YLim',[0 max_velocity])
        hold on
        X_Lims = get(gca,'XLim'); Y_Lims = get(gca,'YLim');
        if L_max == Inf
            L_max = max_length;
        end
        if V_max == Inf
            V_max = max_velocity;
        end
        rectangle('Position',[L_min V_min ...
            L_max-L_min V_max-V_min],...
            'FaceColor','none','EdgeColor',[0.4 0.4 0.4])
        plot(lengths,velocities,'k+')
        if V_thresh > 0
            plot([L_min,L_max],[V_thresh,V_thresh],'k--')
        end
        hold off
        
        % Plot the frame to frame velocity histogram
        subplot(2,1,2)
        
        hist(f2f_velocities,25)
        xlabel('V_{f2f}[\mum/s]')
        ylabel('Count')
        
        if V_thresh > 0
            hold on
            currentLim = get(gca,'YLim');
            plot([V_thresh,V_thresh],[0,currentLim(2).*10.0],'k--')
            set(gca,'YLim',currentLim)
            hold off
        end
        
    end
else
   % If no filament is left over to calculate results for
   
   % Trace velocity quantifiers
    v_mean = NaN; v_STD = NaN;
    L_mean = NaN; L_STD = NaN;
    f_mot = NaN;
   
   fprintf(['\n!!!!!!!\nNo filaments there to calculate results or plot.\n',...
       'Check Trace_Results.dat files in results folders and check\n',...
       'your L_min, L_max, V_min, V_max choices in the interface.\n\n'])
   
end

close(handles.prog_bar)

function write_plain_text(L_min,L_max,V_min,V_max,...
    result_paths,complex_len_in,complex_vel_in,non_dropped,trace_text_file)

% Write filament velocities and lengths to plain text file

number_of_results = numel(result_paths);

lengths = [];
velocities = [];

handles.prog_bar = waitbar(0,'Starting data read in...');

if non_dropped
    % Include all filaments, irrespective of TraceDropper application
    for kk = 1:number_of_results
        
        handles.prog_bar = ...
            waitbar(kk./number_of_results,handles.prog_bar,'Reading in data');
        
        % Extract the trace results for this result file
        this_trace_results = load(result_paths{kk},'trace_results');
        this_trace_results = this_trace_results.trace_results.trace_results;

        this_trace_results = this_trace_results;
        lengths = [lengths,this_trace_results.average_filament_length];
        velocities = [velocities,this_trace_results.trace_velocity];
        
    end
else
    for kk = 1:number_of_results
        % Include only filaments that were found OK in the TraceDropper
        
        handles.prog_bar = ...
            waitbar(kk./number_of_results,handles.prog_bar,'Reading in data');

        
        % Extract the tracedrop-mask for this result file
        tracedrop_mask = load(result_paths{kk},'tracedrop_mask');
        tracedrop_mask = tracedrop_mask.tracedrop_mask;
        
        % Extract the trace results for this result file
        this_trace_results = load(result_paths{kk},'trace_results');
        this_trace_results = this_trace_results.trace_results.trace_results;

        this_trace_results = this_trace_results(tracedrop_mask);
        lengths = [lengths,this_trace_results.average_filament_length];
        velocities = [velocities,this_trace_results.trace_velocity];
        
    end
end

if complex_len_in && complex_vel_in
    in_inds = real(lengths)>=L_min & real(lengths)<=L_max ...
        & real(velocities)>=V_min & real(velocities)<=V_max;
elseif complex_len_in
    in_inds = ~imag(lengths) ...
        & real(lengths)>=L_min & real(lengths)<=L_max ...
        & real(velocities)>=V_min & real(velocities)<=V_max;
elseif complex_vel_in
    in_inds = ~imag(velocities) ...
        & real(lengths)>=L_min & real(lengths)<=L_max ...
        & real(velocities)>=V_min & real(velocities)<=V_max;
else
    in_inds = ~imag(lengths) & ~imag(velocities) ...
        & real(lengths)>=L_min & real(lengths)<=L_max ...
        & real(velocities)>=V_min & real(velocities)<=V_max;
end


% Keep only traces that fulfill the constraints
lengths = lengths(in_inds);
velocities = velocities(in_inds);

% Write length and trace velocity as tab-delimkted text file
field_headers = ...
    {'Length[mum]',...
    'Velocity[mum/s]'};

number_of_columns = numel(field_headers);

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
write_array = cat(2,real(lengths)',real(velocities)');

dlmwrite(trace_text_file,write_array,'-append','delimiter','\t')

close(handles.prog_bar)

function Vmin_edit_Callback(hObject, eventdata, handles)
% hObject    handle to Vmin_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Vmin_edit as text
%        str2double(get(hObject,'String')) returns contents of Vmin_edit as a double


% --- Executes during object creation, after setting all properties.
function Vmin_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Vmin_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Vmax_edit_Callback(hObject, eventdata, handles)
% hObject    handle to Vmax_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Vmax_edit as text
%        str2double(get(hObject,'String')) returns contents of Vmax_edit as a double


% --- Executes during object creation, after setting all properties.
function Vmax_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Vmax_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Lmin_edit_Callback(hObject, eventdata, handles)
% hObject    handle to Lmin_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Lmin_edit as text
%        str2double(get(hObject,'String')) returns contents of Lmin_edit as a double


% --- Executes during object creation, after setting all properties.
function Lmin_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Lmin_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Lmax_edit_Callback(hObject, eventdata, handles)
% hObject    handle to Lmax_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Lmax_edit as text
%        str2double(get(hObject,'String')) returns contents of Lmax_edit as a double


% --- Executes during object creation, after setting all properties.
function Lmax_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Lmax_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Vthresh_edit_Callback(hObject, eventdata, handles)
% hObject    handle to Vthresh_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Vthresh_edit as text
%        str2double(get(hObject,'String')) returns contents of Vthresh_edit as a double


% --- Executes during object creation, after setting all properties.
function Vthresh_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Vthresh_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pickfiles_pushbutton.
function pickfiles_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to pickfiles_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.result_paths = uipickfiles;

set(handles.pickedfiles_text,'String',...
    sprintf('%d result files selected.',numel(handles.result_paths)))

set(handles.files_listbox,'String',handles.result_paths)

guidata(hObject,handles)


% --- Executes on button press in calculate_pushbutton.
function calculate_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to calculate_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.complex_len_in = boolean(get(handles.complex_len_checkbox,'Value'));
handles.complex_vel_in = boolean(get(handles.complex_vel_checkbox,'Value'));
handles.non_dropped = boolean(get(handles.non_dropped_checkbox,'Value'));
handles.do_plots = boolean(get(handles.do_plot_checkbox,'Value'));

handles.V_min = str2double(get(handles.Vmin_edit,'String'));
handles.V_max = str2double(get(handles.Vmax_edit,'String'));

handles.L_min = str2double(get(handles.Lmin_edit,'String'));
handles.L_max = str2double(get(handles.Lmax_edit,'String'));

handles.V_thr = str2double(get(handles.Vthresh_edit,'String'));

[v_mean,v_STD,L_mean,L_STD,f_mot,nn] = ...
    get_motility_quantifiers(...
    handles.L_min,handles.L_max,handles.V_min,handles.V_max,handles.V_thr,...
    handles.result_paths,handles.complex_len_in,handles.complex_vel_in,...
    handles.non_dropped,handles.do_plots);

fprintf('\n\n --- Limits and thresholds --- \n')
fprintf('Length limits: (%f,%f) [microns]\n',...
    handles.L_min,handles.L_max)
fprintf('Velocity limits: (%f,%f) [microns/sec]\n',...
    handles.V_min,handles.V_max)
fprintf('Motile/non-motile threshold velocity: %f [microns/sec]\n',...
    handles.V_thr)

fprintf('\n --- Inclusion options --- \n')
if handles.complex_len_in && handles.complex_vel_in
    fprintf('Complex lengths and velocities included.\n')
elseif handles.complex_len_in
    fprintf('Complex lengths included, complex velocities excluded.\n')
elseif handles.complex_vel_in
    fprintf('Complex velocities included, complex lengths excluded.\n')
else
    fprintf('Complex lengths and velocities excluded.\n')
end

if handles.non_dropped
    fprintf('Trace dropping disregarded.')
else
    fprintf('Dropped traces excluded.')
end

fprintf('\n\n --- Calculated values --- \n')
fprintf('Mean sliding velocity: %f +- %f[microns/sec](StdDev)\n',v_mean,v_STD)
fprintf('Mean filament length: %f +- %f[microns](StdDev)\n',L_mean,L_STD)
fprintf('Motile fraction: %f\n',f_mot)

fprintf('\nFilament count: n = %d, Videos: n = %d\n\n',...
    nn,numel(handles.result_paths))

guidata(hObject,handles)


% --- Executes on button press in complex_len_checkbox.
function complex_len_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to complex_len_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of complex_len_checkbox


% --- Executes on button press in non_dropped_checkbox.
function non_dropped_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to non_dropped_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of non_dropped_checkbox



% --- Executes on button press in do_plot_checkbox.
function do_plot_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to do_plot_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of do_plot_checkbox


% --- Executes on button press in save_pushbutton.
function save_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to save_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.complex_len_in = boolean(get(handles.complex_len_checkbox,'Value'));
handles.complex_vel_in = boolean(get(handles.complex_vel_checkbox,'Value'));
handles.non_dropped = boolean(get(handles.non_dropped_checkbox,'Value'));
handles.do_plots = boolean(get(handles.do_plot_checkbox,'Value'));

handles.V_min = str2double(get(handles.Vmin_edit,'String'));
handles.V_max = str2double(get(handles.Vmax_edit,'String'));

handles.L_min = str2double(get(handles.Lmin_edit,'String'));
handles.L_max = str2double(get(handles.Lmax_edit,'String'));

handles.V_thr = str2double(get(handles.Vthresh_edit,'String'));

[file,path] = uiputfile;
trace_text_file = [path file];

fprintf('Started save to plain text file %s...\n',trace_text_file)

write_plain_text(...
    handles.L_min,handles.L_max,handles.V_min,handles.V_max,...
    handles.result_paths,handles.complex_len_in,handles.complex_vel_in,...
    handles.non_dropped,trace_text_file);

guidata(hObject,handles)


% --- Executes on selection change in files_listbox.
function files_listbox_Callback(hObject, eventdata, handles)
% hObject    handle to files_listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns files_listbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from files_listbox

guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function files_listbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to files_listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in complex_vel_checkbox.
function complex_vel_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to complex_vel_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of complex_vel_checkbox
