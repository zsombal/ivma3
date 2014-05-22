function varargout = AutoDropper(varargin)
% AUTODROPPER MATLAB code for AutoDropper.fig
%      AUTODROPPER, by itself, creates a new AUTODROPPER or raises the existing
%      singleton*.
%
%      H = AUTODROPPER returns the handle to a new AUTODROPPER or the handle to
%      the existing singleton*.
%
%      AUTODROPPER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in AUTODROPPER.M with the given input arguments.
%
%      AUTODROPPER('Property','Value',...) creates a new AUTODROPPER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before AutoDropper_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to AutoDropper_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help AutoDropper

% Last Modified by GUIDE v2.5 21-May-2014 11:57:50

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @AutoDropper_OpeningFcn, ...
                   'gui_OutputFcn',  @AutoDropper_OutputFcn, ...
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


% --- Executes just before AutoDropper is made visible.
function AutoDropper_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to AutoDropper (see VARARGIN)

% Choose default command line output for AutoDropper
handles.output = hObject;

handles.file_path = pwd;
handles.benchmark_flag = false;

handles.threshold = str2double(get(handles.ThresholdEdit,'String'));
handles.threshold = max([0,min([1,handles.threshold])]);
set(handles.ThresholdEdit,'String',sprintf('%3.3f',handles.threshold));
set(handles.ThresholdSlider,'Value',handles.threshold);

disp(handles)

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes AutoDropper wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = AutoDropper_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function TrainingEdit_Callback(hObject, eventdata, handles)
% hObject    handle to TrainingEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TrainingEdit as text
%        str2double(get(hObject,'String')) returns contents of TrainingEdit as a double


% --- Executes during object creation, after setting all properties.
function TrainingEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TrainingEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in TrainingLoadButton.
function TrainingLoadButton_Callback(hObject, eventdata, handles)
% hObject    handle to TrainingLoadButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.TrainingLoadButton,'Enable','off');
set(handles.TrainingCreateButton,'Enable','off');
set(handles.TrainingSaveButton,'Enable','off');

% Which model should be used for tracedropping?
[file_name, handles.file_path] = uigetfile(pwd,...
    'Select tracedrop training file...');
full_model_path = [handles.file_path file_name];

% Load existent model
training_expertise = load(full_model_path);
handles.training_expertise = training_expertise;

handles.decision_tree_model = training_expertise.decision_tree_model;
handles.relevant_pixels = training_expertise.relevant_pixels;
handles.opt_thresh = training_expertise.opt_thresh;
handles.height = training_expertise.height;
handles.width = training_expertise.width;
handles.ntree = handles.decision_tree_model.NTrees;

handles.training_path = full_model_path;

handles.out_of_bag_error = oobError(handles.decision_tree_model);
plot(handles.NumberAxes,1:numel(handles.out_of_bag_error), ...
    handles.out_of_bag_error,'k-')
set(gca,'XLim',[1,numel(handles.out_of_bag_error)])
xlabel('Decision trees in ensemble')
ylabel('Prediction error rate')

set(handles.TrainingLoadButton,'Enable','on');
set(handles.TrainingCreateButton,'Enable','on');
set(handles.TrainingSaveButton,'Enable','on');

set(handles.BenchmarkLoadButton,'Enable','on');
set(handles.BenchmarkSaveButton,'Enable','on');

set(handles.ThresholdSlider,'Enable','on');
set(handles.ThresholdEdit,'Enable','on');
set(handles.ThresholdText,'Enable','on');

set(handles.TrainingEdit,'String',full_model_path);
set(handles.NumberEdit,'String',sprintf('%d',handles.ntree));

set(handles.AutodropButton,'Enable','on');

if isfield(training_expertise,'benchmark_struct')
    benchmark_struct = training_expertise.benchmark_struct;
    handles.false_pos = benchmark_struct.false_pos;
    handles.true_pos = benchmark_struct.true_pos;
    handles.support_scores = benchmark_struct.support_scores;
    handles.threshold = benchmark_struct.threshold;
    handles.benchmark_flag = true;
    update_plots(handles);
    set(handles.BenchmarkEdit,'String','This training file has benchmark');
else
    handles.benchmark_flag = false;
    update_plots(handles);
    set(handles.BenchmarkEdit,'String','No benchmark bundle loaded');
end

set(handles.ImageButton,'Enable','off');

guidata(hObject,handles);


% --- Executes on button press in TrainingCreateButton.
function TrainingCreateButton_Callback(hObject, eventdata, handles)
% hObject    handle to TrainingCreateButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% This code is used to create a training experience file for the automated
% trace dropper system. All parameters for the training have to be
% specified here, they cannot be changed anymore after the training is
% completed.
%
% The only parameter that is still free is the threshold for
% classifying a filament motion trace as rejected, which is between 0 and
% 1. This estimation can be done from the output plots that are presented
% after the training with this file (train_on_bundle.m) has completed.
%
% For background on the classification system, see Hilbert et al., PLoS
% Computation Biology, 2013.

set(handles.TrainingLoadButton,'Enable','off');
set(handles.TrainingCreateButton,'Enable','off');
set(handles.TrainingSaveButton,'Enable','off');

% Algorithm parameters

ntrees = str2num(get(handles.NumberEdit,'String'));

% Open the bundles used for training and testing of the data classification
% system
[file_name, handles.file_path] = uigetfile(handles.file_path,...
    'Select bundle to train classification system.');
training_bundle = [handles.file_path file_name];

non_dropped = true;
sufficient_keywords = {};
necessary_keywords = {};
not_keywords = {};

% ---
% Check if all videos have the same height and width, if not return error
% message

fprintf('Checking if all video dimensions are the same...')

% Get height and width of the different sections
get_height = @(section) ...
    section.video_properties.height;
get_width = @(section) ...
    section.video_properties.width;
heights = ...
    extract_by_keywords( training_bundle, ...
    sufficient_keywords,necessary_keywords,not_keywords, ...
    get_height,non_dropped);
widths = ...
    extract_by_keywords( training_bundle, ...
    sufficient_keywords,necessary_keywords,not_keywords, ...
    get_width,non_dropped);
heights = [heights{:}];
widths = [widths{:}];

fprintf('done.\n')

% Check if heights or widths differ between sections
if numel(unique(heights))>1 || numel(unique(widths))>1 %based on unique elements
    fprintf(['The videos have different dimensions.\n' ...
        'You can not train the system on videos of different dimensions.\n'])
    msgbox(['The videos have different dimensions.' ...
        'You can not train the system on videos of different dimensions.'],...
        'Video dimension mismatch','error')
else
    height = unique(heights);
    width = unique(widths);
end

% ---
% Extract tracedrop masks and trace image data and bounding boxes from
% results bundle

fprintf('Reading in training data...')

get_tracedrop_mask = @(section) ...
    section.tracedrop_mask;
get_trace_image = @(section) ...
    {section.trace_results.trace_image};
get_trace_bounding_box = @(section) ...
    {section.trace_results.trace_bounding_box};
get_trace_length = @(section) ...
    {section.trace_results.trace_length};
get_trace_width = @(section) ...
    {section.trace_results.trace_width};
get_trace_solidity = @(section) ...
    {section.trace_results.trace_solidity};

tracedrop_masks = ...
    extract_by_keywords( training_bundle, ...
    sufficient_keywords,necessary_keywords,not_keywords, ...
    get_tracedrop_mask,non_dropped);
trace_images = ...
    extract_by_keywords( training_bundle, ...
    sufficient_keywords,necessary_keywords,not_keywords, ...
    get_trace_image,non_dropped);
trace_bounding_boxes = ...
    extract_by_keywords( training_bundle, ...
    sufficient_keywords,necessary_keywords,not_keywords, ...
    get_trace_bounding_box,non_dropped);
trace_lengths = ...
    extract_by_keywords( training_bundle, ...
    sufficient_keywords,necessary_keywords,not_keywords, ...
    get_trace_length,non_dropped);
trace_widths = ...
    extract_by_keywords( training_bundle, ...
    sufficient_keywords,necessary_keywords,not_keywords, ...
    get_trace_width,non_dropped);
trace_solidities = ...
    extract_by_keywords( training_bundle, ...
    sufficient_keywords,necessary_keywords,not_keywords, ...
    get_trace_solidity,non_dropped);

tracedrop_masks = [tracedrop_masks{:}];
trace_images = [trace_images{:}];
trace_bounding_boxes = [trace_bounding_boxes{:}];
trace_lengths = [trace_lengths{:}];
trace_lengths = [trace_lengths{:}];
trace_widths = [trace_widths{:}];
trace_widths = [trace_widths{:}];
trace_solidities = [trace_solidities{:}];
trace_solidities = [trace_solidities{:}];

number_of_traces = numel(trace_images);

fprintf('done.\n%d traces in training data set.\n',number_of_traces)

%% Extract features for training data set

fprintf('Extracting trace features from training data set...')

% ---
%Specifications for corner detection
filter_coefficients = fspecial('gaussian',[21 1],2.5);
maximal_corners = 200;

% Determine number of featuresout-of-bundle data from
number_of_features = 2.*maximal_corners + 6;

% Extract and/or store trace features
feature_container = zeros(number_of_traces,number_of_features);
% Create trace images that are centered and have same 
centered_images = false(number_of_traces,height.*width);

fprintf('%d traces to analyze in total. Starting...\n',number_of_traces);

parfor_progress(number_of_traces);

parfor kk = 1:number_of_traces;
    
    % Initialize all 0 image
    this_image = false(height,width);
    
    % horizontal positioning
    xx_center = floor(trace_bounding_boxes{kk}(1) + ...
        trace_bounding_boxes{kk}(3)./2);
    xx_min = ceil(trace_bounding_boxes{kk}(1))-xx_center+round(width./2);
    xx_max = xx_min + trace_bounding_boxes{kk}(3)-1;
    % vertical positioning
    yy_center = floor(trace_bounding_boxes{kk}(2) + ...
        trace_bounding_boxes{kk}(4)./2);
    yy_min = ceil(trace_bounding_boxes{kk}(2))-yy_center+round(height./2);
    yy_max = yy_min + trace_bounding_boxes{kk}(4)-1;
    this_image(yy_min:yy_max,xx_min:xx_max) = ...
        trace_images{kk};
    
    % Get region properties of the trace in this image
    properties = regionprops(this_image,'Area','Orientation','Perimeter');
    
    % Rotate this image to have the trace main axis in horizontal direction
    this_image = imrotate(...
        this_image,-properties.Orientation,'bilinear','crop')
    
    %Store pixel data for this trace
    centered_images(kk,:) = this_image(:); %Store as row in centered_images
    
    %Store feature data for this trace
    corner_positions = corner(this_image,maximal_corners,...
        'FilterCoefficients',filter_coefficients,...
        'SensitivityFactor',0.2,...
        'QualityLevel',0.15);
    detected_corners = numel(corner_positions)./2;
    
    this_trace_features = zeros(1,number_of_features);
    this_trace_features(1:detected_corners) = corner_positions(:,1);
    this_trace_features(...
        (1+maximal_corners):(detected_corners+maximal_corners)) = ...
        corner_positions(:,2);
    % number of pixels in trace
    this_trace_features(2.*maximal_corners+4) = properties.Area;
    % number of detected corners
    this_trace_features(2.*maximal_corners+5) = detected_corners;
    % Perimeter of trace
    this_trace_features(2.*maximal_corners+6) = properties.Perimeter;
    
    feature_container(kk,:) = this_trace_features;
    
    parfor_progress;
    
end
parfor_progress(0);

feature_container(:,2.*maximal_corners+2) = trace_lengths;
feature_container(:,2.*maximal_corners+3) = trace_widths;
feature_container(:,2.*maximal_corners+4) = trace_solidities;

%Remove all imaginary parts from feature values
feature_container = real(feature_container);

% Clear data out of memory
clear trace_images trace_bounding_boxes heights widths

% ---
% Detect pixels for which any variation happens across the training data
% set
fprintf('Determining pixels containing any information...\n')

relevant_pixels = mean(centered_images,1);
relevant_pixels = relevant_pixels > 0 & relevant_pixels<1;

fprintf(...
    '%d of %d pixels contain any information useful for classification...\n', ...
    sum(relevant_pixels),numel(relevant_pixels))

cropped_pixel_data = double(centered_images(:,relevant_pixels));

% Clear data out of memory
clear centered_images

% Merge feature vector and relevant pixels vector
feature_container = cat(2,feature_container,cropped_pixel_data);
clear cropped_pixel_data

fprintf('done.\n')

% ---
% Train bagged decision tree method on training data set

% Create set of aggregated decision trees
fprintf('Creating bootstrap aggregation decision tree model...\n')
opt = statset('UseParallel','always');
decision_tree_model = ...
    TreeBagger(ntrees,feature_container,tracedrop_masks,...
    'oobpred','on','NPrint',1,'Options',opt);
fprintf('done.\n')

% ---
% Make plot showing if enough trees were used
handles.out_of_bag_error = oobError(decision_tree_model);
plot(handles.NumberAxes,1:numel(handles.out_of_bag_error), ...
    handles.out_of_bag_error,'k-')
set(gca,'XLim',[1,numel(handles.out_of_bag_error)])
xlabel('Decision trees in ensemble')
ylabel('Prediction error rate')

% ---
% Store values to handle structure
handles.ntrees = ntrees;
handles.height = height;
handles.width = width;
handles.filter_coefficients = filter_coefficients;
handles.maximal_corners = maximal_corners;
handles.number_of_features = number_of_features;
handles.relevant_pixels = relevant_pixels;
handles.decision_tree_model = decision_tree_model;

handles.training_expertise = struct();
handles.training_expertise.decision_tree_model = handles.decision_tree_model;
handles.training_expertise.relevant_pixels = handles.relevant_pixels;
handles.training_expertise.opt_thresh = 0.8;
handles.training_expertise.height = handles.height;
handles.training_expertise.width = handles.width;

set(handles.TrainingLoadButton,'Enable','on');
set(handles.TrainingCreateButton,'Enable','on');
set(handles.TrainingSaveButton,'Enable','on');

set(handles.ThresholdSlider,'Enable','on');
set(handles.ThresholdEdit,'Enable','on');
set(handles.ThresholdText,'Enable','on');

set(handles.BenchmarkLoadButton,'Enable','on');
set(handles.BenchmarkSaveButton,'Enable','on');

set(handles.TrainingEdit,'String','TRAINING NOT SAVED YET!');

set(handles.AutodropButton,'Enable','on');

set(handles.ImageButton,'Enable','off');

guidata(hObject,handles);


% --- Executes on button press in TrainingSaveButton.
function TrainingSaveButton_Callback(hObject, eventdata, handles)
% hObject    handle to TrainingSaveButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[file_name, handles.file_path] = uiputfile(handles.file_path,...
    'Select file to save trained model to.');
training_file = [handles.file_path, file_name];

decision_tree_model = handles.decision_tree_model;
relevant_pixels = handles.relevant_pixels;
opt_thresh = 0.8;
height = handles.height;
width = handles.width;

set(handles.TrainingLoadButton,'Enable','off');
set(handles.TrainingCreateButton,'Enable','off');
set(handles.TrainingSaveButton,'Enable','off');

save(training_file,'decision_tree_model','relevant_pixels','opt_thresh',...
    'height','width','-v7.3')

handles.training_path = [handles.file_path file_name];
set(handles.TrainingEdit,'String',handles.training_path);

set(handles.TrainingLoadButton,'Enable','on');
set(handles.TrainingCreateButton,'Enable','on');
set(handles.TrainingSaveButton,'Enable','on');

guidata(hObject,handles);

function BenchmarkEdit_Callback(hObject, eventdata, handles)
% hObject    handle to BenchmarkEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of BenchmarkEdit as text
%        str2double(get(hObject,'String')) returns contents of BenchmarkEdit as a double


% --- Executes during object creation, after setting all properties.
function BenchmarkEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to BenchmarkEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in BenchmarkLoadButton.
function BenchmarkLoadButton_Callback(hObject, eventdata, handles)
% hObject    handle to BenchmarkLoadButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[file_name, handles.file_path] = uigetfile(handles.file_path,...
    'Select bundle to benchmark the training file.');
handles.benchmark_path = [handles.file_path file_name];
set(handles.BenchmarkEdit,'String',handles.benchmark_path);

set(handles.BenchmarkLoadButton,'Enable','off');
set(handles.BenchmarkSaveButton,'Enable','off');

% ---
% Test predictive power in second data set


non_dropped = true;
sufficient_keywords = {};
necessary_keywords = {};
not_keywords = {};

% ---
% Check if all videos have the same height and width, if not return error
% message

fprintf('Checking if all video dimensions are the same...')

% Get height and width of the different sections
get_height = @(section) ...
    section.video_properties.height;
get_width = @(section) ...
    section.video_properties.width;
heights = ...
    extract_by_keywords( handles.benchmark_path, ...
    sufficient_keywords,necessary_keywords,not_keywords, ...
    get_height,non_dropped);
widths = ...
    extract_by_keywords( handles.benchmark_path, ...
    sufficient_keywords,necessary_keywords,not_keywords, ...
    get_width,non_dropped);
heights = [heights{:}];
widths = [widths{:}];

fprintf('done.\n')

% Check if heights or widths differ between sections
disp(handles)
if numel(unique(heights))>1 || numel(unique(widths))>1 ...
        || unique(heights) ~= handles.height ...
        || unique(widths) ~= handles.width %based on unique elements
    fprintf(['The videos have different dimensions.\n' ...
        'You can not benchmark the system on videos of different dimensions.\n'])
    msgbox(['The videos have different dimensions.' ...
        'You can not benchmark the system on videos of different dimensions.'],...
        'Video dimension mismatch','error')
end


% Extract tracedrop masks and trace image data and bounding boxes from
% results bundle

fprintf('Reading in test data...')

get_tracedrop_mask = @(section) ...
    section.tracedrop_mask;
get_trace_image = @(section) ...
    {section.trace_results.trace_image};
get_trace_bounding_box = @(section) ...
    {section.trace_results.trace_bounding_box};
get_trace_length = @(section) ...
    {section.trace_results.trace_length};
get_trace_width = @(section) ...
    {section.trace_results.trace_width};
get_trace_solidity = @(section) ...
    {section.trace_results.trace_solidity};

tracedrop_masks = ...
    extract_by_keywords( handles.benchmark_path, ...
    sufficient_keywords,necessary_keywords,not_keywords, ...
    get_tracedrop_mask,non_dropped);
trace_images = ...
    extract_by_keywords( handles.benchmark_path, ...
    sufficient_keywords,necessary_keywords,not_keywords, ...
    get_trace_image,non_dropped);
trace_bounding_boxes = ...
    extract_by_keywords( handles.benchmark_path, ...
    sufficient_keywords,necessary_keywords,not_keywords, ...
    get_trace_bounding_box,non_dropped);
trace_lengths = ...
    extract_by_keywords( handles.benchmark_path, ...
    sufficient_keywords,necessary_keywords,not_keywords, ...
    get_trace_length,non_dropped);
trace_widths = ...
    extract_by_keywords( handles.benchmark_path, ...
    sufficient_keywords,necessary_keywords,not_keywords, ...
    get_trace_width,non_dropped);
trace_solidities = ...
    extract_by_keywords( handles.benchmark_path, ...
    sufficient_keywords,necessary_keywords,not_keywords, ...
    get_trace_solidity,non_dropped);

tracedrop_masks = [tracedrop_masks{:}];
trace_images = [trace_images{:}];
trace_bounding_boxes = [trace_bounding_boxes{:}];
trace_lengths = [trace_lengths{:}];
trace_lengths = [trace_lengths{:}];
trace_widths = [trace_widths{:}];
trace_widths = [trace_widths{:}];
trace_solidities = [trace_solidities{:}];
trace_solidities = [trace_solidities{:}];

number_of_traces = numel(trace_images);

fprintf('Extracting trace features from test data set...')

% ---
%Specifications for corner detection
filter_coefficients = fspecial('gaussian',[21 1],2.5);
maximal_corners = 200;

% Determine number of featuresout-of-bundle data from
number_of_features = 2.*maximal_corners + 6;

% Extract and/or store trace features
feature_container = zeros(number_of_traces,number_of_features);
% Create trace images that are centered and have same 
centered_images = false(number_of_traces,handles.height.*handles.width);

parfor_progress(number_of_traces);

parfor kk = 1:number_of_traces;
    
    % Initialize all 0 image
    this_image = false(handles.height,handles.width);
    
    % horizontal positioning
    xx_center = floor(trace_bounding_boxes{kk}(1) + ...
        trace_bounding_boxes{kk}(3)./2);
    xx_min = ceil(trace_bounding_boxes{kk}(1))-xx_center+round(handles.width./2);
    xx_max = xx_min + trace_bounding_boxes{kk}(3)-1;
    % vertical positioning
    yy_center = floor(trace_bounding_boxes{kk}(2) + ...
        trace_bounding_boxes{kk}(4)./2);
    yy_min = ceil(trace_bounding_boxes{kk}(2))-yy_center+round(handles.height./2);
    yy_max = yy_min + trace_bounding_boxes{kk}(4)-1;
    this_image(yy_min:yy_max,xx_min:xx_max) = ...
        trace_images{kk};
    
    % Get region properties of the trace in this image
    properties = regionprops(this_image,'Area','Orientation','Perimeter');
    
    % Rotate this image to have the trace main axis in horizontal direction
    this_image = imrotate(...
        this_image,-properties.Orientation,'bilinear','crop')
    
    %Store pixel data for this trace
    centered_images(kk,:) = this_image(:); %Store as row in centered_images
    
    %Store feature data for this trace
    corner_positions = corner(this_image,maximal_corners,...
        'FilterCoefficients',filter_coefficients,...
        'SensitivityFactor',0.2,...
        'QualityLevel',0.15);
    detected_corners = numel(corner_positions)./2;
    
    this_trace_features = zeros(1,number_of_features);
    this_trace_features(1:detected_corners) = corner_positions(:,1);
    this_trace_features(...
        (1+maximal_corners):(detected_corners+maximal_corners)) = ...
        corner_positions(:,2);
    % number of pixels in trace
    this_trace_features(2.*maximal_corners+4) = properties.Area;
    % number of detected corners
    this_trace_features(2.*maximal_corners+5) = detected_corners;
    % Perimeter of trace
    this_trace_features(2.*maximal_corners+6) = properties.Perimeter;
    
    feature_container(kk,:) = this_trace_features;
    
    parfor_progress;
    
end

parfor_progress(0);

feature_container(:,2.*maximal_corners+2) = trace_lengths;
feature_container(:,2.*maximal_corners+3) = trace_widths;
feature_container(:,2.*maximal_corners+4) = trace_solidities;

%Remove all imaginary parts from feature values
feature_container = real(feature_container);

% Clear data out of memory
clear trace_images trace_bounding_boxes heights widths

cropped_pixel_data = double(centered_images(:,handles.relevant_pixels));

% Clear data out of memory
clear centered_images

% Merge feature vector and relevant pixels vector
feature_container = cat(2,feature_container,cropped_pixel_data);
clear cropped_pixel_data

fprintf('done.\n')


% ---
% Make class predictions in test data set
fprintf('Making class predictions for test data set...')
%plot(oobError(decision_tree_model))
[predicted_class,prediction_score] = ...
    predict(handles.decision_tree_model,feature_container);
predicted_class = logical(cellfun(@(xx) str2num(xx),predicted_class));
score_value = prediction_score(:,2);
fprintf('done.\n')


% True vs. false positive ROC curve
[false_pos,true_pos,support_scores] = ...
    perfcurve(tracedrop_masks,score_value,true);

axes(handles.ROCAxes);
plot(handles.ROCAxes,false_pos,true_pos,'k-','LineWidth',0.75)
xlabel('False positive rate')
ylabel('True positive rate')
legend('ROC curve','location','southeast')
legend BOXOFF

axes(handles.MisClassAxes);
plot(handles.MisClassAxes,support_scores,true_pos,'k-','LineWidth',0.75)
hold on
plot(handles.MisClassAxes,support_scores,1.0-false_pos,'b--','LineWidth',0.75)
hold off
xlabel('Threshold')
ylabel('Correct fraction')
legend('Keep','Reject','location','Northwest')
legend BOXOFF;
set(gca,'YLim',[0,1],'XLim',[0,1])

handles.score_value = score_value;
handles.false_pos = false_pos;
handles.true_pos = true_pos;
handles.support_scores = support_scores;
handles.tracedrop_masks = tracedrop_masks;

handles.benchmark_flag = true;

set(handles.BenchmarkLoadButton,'Enable','on');
set(handles.BenchmarkSaveButton,'Enable','on');

set(handles.AutodropButton,'Enable','on');

set(handles.ImageButton,'Enable','on');
set(handles.ImageText,'Enable','on');

update_plots(handles)

guidata(hObject,handles);

% --- Executes on button press in BenchmarkSaveButton.
function BenchmarkSaveButton_Callback(hObject, eventdata, handles)
% hObject    handle to BenchmarkSaveButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Add a structure to the saved training file

benchmark_struct = struct();
benchmark_struct.false_pos = handles.false_pos;
benchmark_struct.true_pos = handles.true_pos;
benchmark_struct.support_scores = handles.support_scores;
benchmark_struct.threshold = handles.threshold;

save(handles.training_path,'benchmark_struct','-append');

guidata(hObject,handles);

% --- Executes on slider movement.
function ThresholdSlider_Callback(hObject, eventdata, handles)
% hObject    handle to ThresholdSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

handles.threshold = get(hObject,'Value');
set(handles.ThresholdEdit,'String',sprintf('%3.3f',handles.threshold));

update_plots(handles);

guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function ThresholdSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ThresholdSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function ThresholdEdit_Callback(hObject, eventdata, handles)
% hObject    handle to ThresholdEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ThresholdEdit as text
%        str2double(get(hObject,'String')) returns contents of ThresholdEdit as a double

handles.threshold = str2double(get(hObject,'String'));
handles.threshold = max([0,min([1,handles.threshold])]);
set(hObject,'String',sprintf('%3.3f',handles.threshold));
set(handles.ThresholdSlider,'Value',handles.threshold);

update_plots(handles);

guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function ThresholdEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ThresholdEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function TraceSlider_Callback(hObject, eventdata, handles)
% hObject    handle to TraceSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


slider_val = get(hObject,'Value');
slider_val = round(slider_val);
set(hObject,'Value',slider_val);

% Initialize all 0 image
this_image = false(handles.height,handles.width);

kk = slider_val;

% horizontal positioning
xx_center = floor(handles.trace_bounding_boxes{kk}(1) + ...
    handles.trace_bounding_boxes{kk}(3)./2);
xx_min = ceil(handles.trace_bounding_boxes{kk}(1))-xx_center+round(handles.width./2);
xx_max = xx_min + handles.trace_bounding_boxes{kk}(3)-1;
% vertical positioning
yy_center = floor(handles.trace_bounding_boxes{kk}(2) + ...
    handles.trace_bounding_boxes{kk}(4)./2);
yy_min = ceil(handles.trace_bounding_boxes{kk}(2))-yy_center+round(handles.height./2);
yy_max = yy_min + handles.trace_bounding_boxes{kk}(4)-1;
this_image(yy_min:yy_max,xx_min:xx_max) = ...
    handles.trace_images{kk};

colormap(gray)
axes(handles.TraceAxes)
imagesc(this_image)
    
guidata(hOBject,handles);


% --- Executes during object creation, after setting all properties.
function TraceSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TraceSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in ImageButton.
function ImageButton_Callback(hObject, eventdata, handles)
% hObject    handle to ImageButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

sufficient_keywords = {};
necessary_keywords = {};
not_keywords = {};
non_dropped = true;


fprintf('Loading missed traces'' images...')
get_trace_image = @(section) ...
    {section.trace_results.trace_image};
get_trace_bounding_box = @(section) ...
    {section.trace_results.trace_bounding_box};
trace_images = ...
    extract_by_keywords( handles.benchmark_path, ...
    sufficient_keywords,necessary_keywords,not_keywords, ...
    get_trace_image,non_dropped);
handles.trace_images = [trace_images{:}];
trace_bounding_boxes = ...
    extract_by_keywords( handles.benchmark_path, ...
    sufficient_keywords,necessary_keywords,not_keywords, ...
    get_trace_bounding_box,non_dropped);
handles.trace_bounding_boxes = [trace_bounding_boxes{:}];
fprintf('done.\n')

handles.false_positive_indices = find( ...
    handles.score_value>=handles.threshold ...
    & handles.tracedrop_masks.'==false);

set(handles.TraceSlider,'Min',1, ...
    'Max',numel(handles.false_positive_indices), ...
    'Value',1,'Enable','on', ...
    'SliderStep',[1./numel(handles.false_positive_indices), ...
    10./numel(handles.false_positive_indices)]);

guidata(hObject,handles);

function NumberEdit_Callback(hObject, eventdata, handles)
% hObject    handle to NumberEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of NumberEdit as text
%        str2double(get(hObject,'String')) returns contents of NumberEdit as a double


% --- Executes during object creation, after setting all properties.
function NumberEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NumberEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function update_plots(handles)


if handles.benchmark_flag
    
    axes(handles.ROCAxes);
    plot(handles.ROCAxes,handles.false_pos,handles.true_pos, ...
        'k-','LineWidth',0.75)
    hold on
    [unique_scores,unique_inds,~] = unique(handles.support_scores);
    unique_true_pos = handles.true_pos(unique_inds);
    unique_false_pos = handles.false_pos(unique_inds);
    true_pos_val = ...
        interp1(unique_scores,unique_true_pos,handles.threshold);
    false_pos_val = ...
        interp1(unique_scores,unique_false_pos,handles.threshold);
    plot([0,1],[1,1]*true_pos_val,'r--')
    plot([1,1]*false_pos_val,[0,1],'r--')
    hold off
    xlabel('False positive rate')
    ylabel('True positive rate')
    legend('ROC curve','location','southeast')
    legend BOXOFF
    
    axes(handles.MisClassAxes);
    plot(handles.MisClassAxes,handles.support_scores,handles.true_pos, ...
        'k-','LineWidth',0.75)
    hold on
    plot(handles.MisClassAxes,handles.support_scores,1.0-handles.false_pos, ...
        'b--','LineWidth',0.75)
    plot([1,1]*handles.threshold,[0,1],'r--')
    hold off
    xlabel('Threshold')
    ylabel('Correct fraction')
    legend('Keep','Reject','location','Northwest')
    legend BOXOFF;
    set(gca,'YLim',[0,1],'XLim',[0,1])
    
    set(handles.KeepText,'String', ...
        sprintf('Keep: %1.1f %% correct',100.*true_pos_val));
    set(handles.RejectText,'String', ...
        sprintf('Reject: %1.1f %% correct',100.*(1.0-false_pos_val)));
    
end


% --- Executes on button press in AutodropButton.
function AutodropButton_Callback(hObject, eventdata, handles)
% hObject    handle to AutodropButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Which directory should be tracedropped?
parent_directory = uigetdir(handles.file_path,...
    'Choose root directory for automated trace dropping.');

% Overwrite existent tracedrop_masks?
overwrite = 'Yes';

% ---
% Call function that automatically drops traces for all
% Analysis_Results.mat files in parent_directory and its subfolders

% Use optimal threshold in detection score that is saved with the model
% elapsed_time = ...
%     auto_tracedrop(training_expertise,parent_directory,overwrite);

% Use user-specified threshold for detection score
elapsed_time = ...
    auto_tracedrop(handles.training_expertise,parent_directory, ...
    overwrite,handles.threshold);

fprintf('Elapsed time: %f seconds.\n',elapsed_time);

guidata(hObject,handles);
