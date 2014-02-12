function [ analysis_output,window_centers,window_counts,varargout ] = ...
    apply_length_resolved( ...
    input_bundle,keywords,L_min,L_max,windows,window_width,...
    complex_in,non_dropped, ...
    filament_function,window_function,varargin)
%APPLY_LENGTH_RESOLVED Applies analysis function in length-resolved manner
%
%   This is a utility function to interact with result bundles. The
%   filaments will first be analyzed individually, then they are binned
%   into running filament length windows, and a summary statistic or result
%   is calculated in each length window.
%
%   The general idea is:
%
%     (1) Specify the result bundle(s) that should be analyzed
%     (2) Specify the function(s) that should be applied for the analysis
%     of each individual filament
%     (3) Specify the function(s) that will be applied to the results of
%     the filaments contained within a given length window.
%
%   For several inputs and outputs, they can be passed in different
%   formats. Also, some inputs and outputs are optional, meaning that they
%   can be used or omitted dependent on what type of analysis is required.
%
%   INPUTS
%
%   input_bundle - bundle from which result sections are read, either
%   passed as a string, or as a cell of strings
%
%   keywords - struct_array with fields and_keyword, with_keywords,
%   not_keywords. If empty [], no keywords are applied
%
%   L_min - start point for box car run
%   L_max - end point for box car run
%   windows - number of box car windows
%   window_width - width of box car window in fraction of the distance from
%   L_min to L_max
%
%   complex_in - if true, complex results in length and velocity are
%   considered, otherwise they are not considered
%
%   non_dropped - is assigned true, all filaments are considered
%   irrespective of if they were removed using the Tracedropper or not. If
%   anything other than true, i.e. unassigned, empty [], or false, the
%   Tracedropper rejection is considered
%
%   filament_function - function to get the necessary results from each
%   filament
%   window_function - function to be applied to all results inside each
%   window
%
%   OPTIONAL INPUTS - if not used, either omit, or use an empty array [] as
%   the argument. The optional arguments have to be added in the right
%   order. So if the third optional argument should be passed, the first
%   and second have to be passed as well. However, just assing the empty
%   array [] in their place if they are not to be applied.
%
%   verbose (first optional argument) - if the string 'verbose' is used as
%   the first optional input argument, the function will output additional
%   details on its progress. Otherwise minimal output will be generated.
%
%   selection_function (second optional argument) - function handle to a 
%   function that is applied to an individual section of the result bundle.
%   If the function returns true, the result section will be included in
%   the analysis, otherwise it will be excluded.
%
%   nboot - number of bootstrap resamples. Only relevant if bootstrap
%   sampling has been requested by using optional output arguments.
%
%   OUTPUTS
%
%   analysis_output - returns the length-resolved result from the window
%   function. The format returned is a 1-by-M cell, where M is the number of
%   level functions, with each element containing a 1-by-N cell, where N is
%   the number of length windows.
%
%   window_centers - returns the centers of the filament length windows,
%   the format returned is a 1-by-N array, where N is the number of length
%   windows.
%
%   window_coutns - returns the number of filaments contained within each
%   length window.
%
%   OPTIONAL OUTPUTS - if not used, either omit, or use tilde ~ as
%   the output. The optional outputs have to be added in the right
%   order. So if the third optional output should be assigned, the first
%   and second have to be passed as well. However, just assign tilde ~
%   in their place if they are not to be assigned.
%
%   confidence_intervals - returns the confidence intervals of the
%   statistic, the confidence intervals are determined from bootstrap
%   statistics, executed with nboot bootstrap samples (500 default value
%   for nboot, if not chosen via optional input argument). The confidence
%   intervals are returned in a cell array of dimension 1-by-M (where M is
%   the number of window level analysis functions). Each cell element
%   contains a 2-by-N numeric array, where the first dimension refers to
%   the upper and lower confidence interval, and the second dimension
%   refers to the actin length window.
%
%   bootstrap_results - returns the bootstrap results that were used in the
%   determination of the confidence intervals. The bootstrap samples
%   are returned in a cell array of dimension 1-by-M (where M is
%   the number of window level analysis functions). Each cell element
%   contains an nbbot-by-N numeric array, where the first dimension refers to
%   the individual bootstrap samples, and the second dimension
%   refers to the actin length window. For one index in the first
%   dimension, first bootstrapping is applied, and then the running window
%   procedure is applied. In this way, for the same value for the first
%   dimension index, all values for the length windows stem from the same
%   re-sample.
%
%   Confidence intervals only work where the window level analysis
%   functions return scalar values.


cell(1,number_of_win_functions);
       for ff = 1:number_of_win_functions
           starting_ind = 1 + (ff-1).*windows;
           varargout{1}{ff} = ...
               boot_ci(:,starting_ind:starting_ind+windows-1)



%% Check if verbose progress report is required

if ~isempty(varargin) && strcmp(varargin{1},'verbose')
    verbose = true;
else
    verbose = false;
end

%% Check if a selection function should be applied to eeach section

if nargin > 11 && ~isempty(varargin{2}) ...
        && strcmp(class(varargin{2}),'function_handle')
    % Apply a selection function to each section
    selection_function = varargin{2};
else
    selection_function = []; % Empty selection function leads to no selection
end



%% Parameters

if nargin > 12 && ~isempty(varargin{3}) ...
        && isnumeric(varargin{3})
    % Set the number of bootstrap samples according to function argument
    nboot = round(varargin{3});
else
    nboot = 500; % Default number of bootstrap samples
end



%% Preliminaries

% ---
% Input bundle acquisition

if verbose
    fprintf('Acquiring result bundles...\n')
end

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
    is_mat_file = exist(input_bundle,'file');
    
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
    L_max = 14;
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
if ~exist('window_width','var') || isempty(window_width) || ...
        (window_width <= 0 || window_width >= 1)
    window_width = 0.2;
end


if ~exist('non_dropped','var') || isempty(non_dropped) || ...
        ~non_dropped
    non_dropped = false;
end

%% Check that functions were passed correctly

if ~iscell(filament_function)
    filament_function = {filament_function};
    filfun_iscell = false;
else
    filfun_iscell = true;
end
number_of_fil_functions = numel(filament_function);
if ~iscell(window_function)
    window_function = {window_function};
    winfun_iscell = false;
else
    winfun_iscell = true;
end
number_of_win_functions = numel(window_function);
if winfun_iscell ~= filfun_iscell
   error(['Filament and window level functions have been passed as' ...
       'different array types. Make both cell array or plain function handle.']) 
end
if ~winfun_iscell && ~filfun_iscell
    % If the filament as well as the window level function are both passed
    % as plain function handles, return the results not as a cell array,
    % but as a plain numeric array
    return_plain = true;
end
if number_of_fil_functions ~= number_of_win_functions
   error('Different number of filament and window level functions.') 
end



%% Get average filament lengths, trace velocities, filament level results

if verbose
    fprintf('Extracting filament velocities...\n')
end

% Pull trace velocities out of the queried result sections
pull_velocity = ...
    @(result_section) [result_section.trace_results.trace_velocity];
trace_velocity = extract_by_keywords(input_bundle, ...
    keywords.any_keywords,keywords.with_keywords,keywords.not_keywords, ...
    pull_velocity,non_dropped,selection_function);

if verbose
    fprintf('Extracting filament lengths...\n')
end

% Pull filament average lenghts out of the queried result sections
pull_average_length = ...
    @(result_section) ...
    [result_section.trace_results.average_filament_length];
average_filament_length = extract_by_keywords(input_bundle, ...
    keywords.any_keywords,keywords.with_keywords,keywords.not_keywords, ...
    pull_average_length,non_dropped,selection_function);

if verbose
    fprintf('Extracting filament results...\n')
end

% Apply the filament level analysis functions to get filament level results
filament_results = cell(1,number_of_fil_functions);
for ff = 1:number_of_fil_functions
    filament_results{ff} = extract_by_keywords(input_bundle, ...
        keywords.any_keywords,keywords.with_keywords,keywords.not_keywords, ...
        filament_function{ff},non_dropped,selection_function);
end

% Remove complex solutions if requested
number_of_sections = numel(average_filament_length);
for cc = 1:number_of_sections
    if ~complex_in
       inds = find(imag(average_filament_length{cc})==0 ...
           & imag(trace_velocity{cc})==0);
       average_filament_length{cc} = average_filament_length{cc}(inds);
       trace_velocity{cc} = trace_velocity{cc}(inds);
       for ff = 1:number_of_fil_functions
           filament_results{ff}{cc} = filament_results{ff}{cc}(inds);
       end
    end
end

% Pool results across all sections
pooled_trace_velocity = [trace_velocity{:}];
pooled_average_filament_length = [average_filament_length{:}];
pooled_filament_results = ...
    cellfun(@(results)[results{:}],filament_results,...
    'UniformOutput',false);

%% Execute the box averaging procedure

if verbose
    fprintf('Calculating length window level features...\n')
end

% ---
% Specify the window lower and upper bounds

interval_length = L_max-L_min; % Length of whole interval
absolute_window_width = interval_length.*window_width;
window_lower_bounds = ...
    linspace(L_min,L_max-absolute_window_width,windows);
window_upper_bounds = ...
    linspace(L_min+absolute_window_width,L_max,windows);
window_centers = ...
    window_lower_bounds + absolute_window_width./2;

% ---
% Execute window level function to calculate window results

% Allocate array to contain the analysis results
analysis_output = cell(1,number_of_win_functions);

% Allocate array to count the filaments in each window
window_counts = zeros(1,windows);
for ff = 1:number_of_win_functions
    analysis_output{ff} = cell(1,windows);
end

for ww = 1:windows
   window_inds = ...
       pooled_average_filament_length>=window_lower_bounds(ww) & ...
       pooled_average_filament_length<=window_upper_bounds(ww);
   for ff = 1:number_of_win_functions
       analysis_output{ff}{ww} = ...
           window_function{ff}(pooled_filament_results{ff}(window_inds));
   end
   window_counts(ww) = numel(window_inds);
end

if return_plain
    analysis_output = analysis_output{1};
end

%% Calculation of bootstrap confidence intervals for window function


if nargout > 3
   % If confidence intervals are requested as a function output, bootstrap
   % results by including different sections

   if verbose
       fprintf('Determining bootstrap confidence intervals...\n')
   end

   boot_function = @(indices) ...
       boot_boxcar_analysis(windows,window_lower_bounds, ...
       window_upper_bounds,window_function,...
       average_filament_length,filament_results,indices);
      
   if nargout > 4
       varargout = cell(1,2);
       [boot_ci,bootstats] = bootci(nboot,boot_function,...
           1:number_of_sections);
       varargout{1} = cell(1,number_of_win_functions);
       varargout{2} = cell(1,number_of_win_functions);
       for ff = 1:number_of_win_functions
           starting_ind = 1 + (ff-1).*windows;
           varargout{1}{ff} = ...
               boot_ci(:,starting_ind:starting_ind+windows-1);
           varargout{2}{ff} = ...
               bootstats(:,starting_ind:starting_ind+windows-1);
       end
       if return_plain
           varargout{1} = varargout{1}{1};
           varargout{2} = varargout{2}{1};
       end

   else
       [boot_ci] = bootci(nboot,boot_function,...
           1:number_of_sections);
       varargout{1} = cell(1,number_of_win_functions);
       for ff = 1:number_of_win_functions
           starting_ind = 1 + (ff-1).*windows;
           varargout{1}{ff} = ...
               boot_ci(:,starting_ind:starting_ind+windows-1);
       end
       if return_plain
           varargout{1} = varargout{1}{1};
       end

   end
   
end


function analysis_output = ...
    boot_boxcar_analysis(windows,window_lower_bounds, ...
    window_upper_bounds,window_function,...
    filament_length,filament_results,indices)


number_of_win_functions = numel(window_function);


% Pool data from sections according to indices passed by bootstrap
% procedure
filament_length = [filament_length{indices}];
for ff = 1:number_of_win_functions
    filament_results{ff} = [filament_results{ff}{indices}];
end


% Allocate array to contain the analysis results
analysis_output = zeros(windows,number_of_win_functions);
% Generate combined feature extraction function
window_function_implicit = @(filament_results,window_inds) ...
    combined_window_function(window_function,filament_results,window_inds);

for ww = 1:windows
    window_inds = ...
        filament_length>=window_lower_bounds(ww) & ...
        filament_length<=window_upper_bounds(ww);
    analysis_output(ww,:) = ...
        window_function_implicit(filament_results,window_inds);
end


function window_results = ...
    combined_window_function(window_function,filament_results,window_inds)

number_of_win_functions = numel(window_function);
window_results = zeros(1,number_of_win_functions);
for ff = 1:number_of_win_functions
    window_results(ff) = ...
        window_function{ff}(filament_results{ff}(window_inds));
end