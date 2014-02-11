function [ output_args ] = mean_trace_velocity( ...
    input_bundle,keywords,L_min,L_max,windows,window_width,,complex_in)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%   Inputs
%
%   input_bundle - bundle from which result sections are read
%
%   keywords - struct_array with fields and_keyword, with_keywords,
%   not_keywords. If empty [], no keywords are applied
%
%   complex_in - if true, complex results in length and velocity are
%   considered, otherwise they are not considered
%
%   L_min - start point for box car run
%   L_max - end point for box car run
%   windows - number of box car windows
%   window_width - width of box car window in fraction of the distance from
%   L_min to L_max



%% Preliminaries

% ---
% Input bundle acquisition

if ~exist('input_bundle','var')
    input_bundle = [];
end

if isempty(input_bundle)
    % If an empty array is given as input_handle, choose an input bundle
    % from disk
    [load_file,load_path] = ...
        uigetfile(pwd,'Select result bundle file from disk.','*.mat');
    input_bundle = [load_path load_file];
    
elseif ischar(input_bundle)
    % If a string is given as the input bundle, check if the string points
    % towards a directory, and if so let user pick a result bundle file
    % starting from that directory. If the string points towards a
    % .mat-file, try to use this file as a result bundle file to load
    
    if exist(input_bundle,'file') && is_mat_file==7 %means that path is a directory
        % Pick file dialog with pre-set start directory
        [load_file,load_path] = ...
            uigetfile(input_bundle,'Select result bundle file from disk.');
        input_bundle = [load_path load_file];
    end
        
end

% ---
% Make sure that keyword fields are functional

if ~exist('keywords','var')
    keywords = [];
end

if ~isstruct(keywords)
    keywords = struct;
    keywords.any_keywords = [];
    keywords.with_keywords = [];
    keywords.not_keywords = [];
else
    if ~isfield(keywords,'and_keywords')
        keywords.and_keywords = [];
    end
    if ~isfield(keywords,'with_keywords')
        keywords.with_keywords = [];
    end
    if ~isfield(keywords,'not_keywords')
        keywords.not_keywords = [];
    end
end

% ---
% Assign complex_in with false if it is anything but true
if ~exist('complex_in','var') || isempty(complex_in) || ~complex_in
    complex_in = false;
else
    complex_in = true;
end

% ---
% Assign parmaeters if they are not given as function inputs, and make sure
% they have functional values
if ~exist('L_min','var') || isempty(L_min) 
    L_min = 0;
end
if ~exist('L_max','var') || isempty(L_max)
    L_max = Inf;
end
if L_min == L_max
    %if L_min is the same as L_max, the window has to be widened
    fprintf('L_min and L_max are the same, L_min lowered to 0.\n')
    L_min = 0;
end
if L_min > L_max
    % Inverse L_min and L_max so that L_min is smaller than L_max
    fprintf('L_max <L_min, reversed to fix it.\n')
    [L_max,L_min] = deal(L_min,L_max);
end
if ~exist('windows','var') || isempty(windows) || windows < 1
    windows = 8;
end
if ~exist('windows','var') || isempty(windows) || windows < 1
    windows = 0.2;
end




%% Get the average filament lengths and the trace velocities

% Pull trace velocities out of the queried result sections
pull_velocity = ...
    @(result_section) [result_section.trace_results( ...
    result_section.tracedrop_mask...
    ).trace_velocity];
trace_velocity = extract_by_keywords(input_bundle, ...
    keywords.any_keywords,keywords.with_keywords,keywords.not_keywords, ...
    pull_velocity);
trace_velocity = [trace_velocity{:}];

% Pull filament average lenghts out of the queried result sections
pull_average_length = ...
    @(result_section) ...
    [result_section.trace_results( ...
    result_section.tracedrop_mask...
    ).average_filament_length];
average_filament_length = extract_by_keywords(input_bundle, ...
    keywords.any_keywords,keywords.with_keywords,keywords.not_keywords, ...
    pull_average_length);
average_filament_length = [average_filament_length{:}];


%% Apply filament length limits and other selection criteria

% Apply the L_min and L_max limits
in_limits_inds = find( ...
    average_filament_length >= L_min ...
    & average_filament_length <= L_max);

% Find indices for which the filament length has no imaginary part, i.e.
% the indices of all filaments for which the rectangular transformation has
% a real solution
non_imag_inds = in_limits_inds(...
    imag(average_filament_length(in_limits_inds))==0 ...
    & imag(trace_velocity)==0);

if complex_in
    inds = in_limits_inds;
else
    inds = non_imag_inds
end

% Plot the (V,L) diagram
VL_figure = figure;
VL_axes = axes('Parent',VL_figure);
plot(VL_axes, ...
    real(average_filament_length(inds)), ...
    real(trace_velocity(inds)), ...
    'ko','MarkerSize',4)
set(VL_axes,'NextPlot','Add')
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
