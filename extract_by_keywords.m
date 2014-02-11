function extracted_results = ...
    extract_by_keywords( input_bundle, ...
    sufficient_keywords,necessary_keywords,not_keywords, ...
    query_function,non_dropped,varargin)

%EXTRACT_BY_KEYWORDS Utility function to extract features from one or
%several result bundles from motility analysis, filter keywords can be
%applied
%
% extracted_results = extract_by_keywords( ...
%   sufficient_keywords,necessary_keywords,not_keywords, ...
%   query_function,non_dropped)
%
%   input_bundle
%   is a cell array containing the structure arrays of the
%   bundled result sections. If an empty array is passed as the bundle, the
%   user is presented with a file selection dialog allowing to
%   interactively select a bundle file to be loaded, if a string is passed,
%   the string will be used to either provide a start path for the file
%   selection dialog, or to directly load a result bundle file if the
%   string points towards a .mat file
%
%   sufficient_keywords
%   is a cell array of strings containing all keywords
%   that by themselves select a result section. Pass [] for no keywords.
%
%   necessary_keywords
%   is a cell array of strings containing all keywords
%   which a section must have to be extracted. Pass [] for no keywords.
%
%   not_keywords
%   is a cell array of strings that must not occur in the
%   keyword list. Pass [] for no keywords.
%
%   non_dropped
%   if assigned true, all filaments are considered
%   irrespective of if they were removed using the Tracedropper or not. If
%   anything other than true, i.e. unassigned, empty [], or false, the
%   Tracedropper rejection is considered
%
%   extracted_results
%   extracted_results is a cell array
%   containing the result of the query_function applied to each section
%   that matches the keyword requirements.
%
%   Example:
%   
%   % Code to plot a histogram of filament lengths in a result bundle
%   non_dropped = false;
%   % Choose an input bundle from disk
%   [load_file,load_path] = ...
%       uigetfile('Select result bundle file from disk.');
%   result_bundle = [load_path load_file];
%   % Define a function that extracts a specific feature from each section
%   % of the result bundle. Here, all filaments' lengths are extracted
%   get_length = @(section) ...
%       [section.trace_results.average_filament_length];
%   % Extract feature from each section, without any filter keywords, note
%   % that the results arrive in a cell array, one cell per section in the 
%   % result bundle
%   lengths = extract_by_keywords(result_bundle,...
%       [],[],[],get_length,non_dropped);
%   % Here the lengths from the cell elements are merged into one regular 
%   % numeric array
%   lengths = [lengths{:}];
%   hist(lengths)
%
%
%  ------------------
%  Optional argument:
%
%  extracted_results = extract_by_keywords( ...
%    sufficient_keywords,necessary_keywords,not_keywords, ...
%    query_function,non_dropped,selection_function)
%
%  selection-function
%  function handle that contains a function to determine if the bundle
%  fulfills a certain selection criterion. The selection function must have
%  a single input argument, which is a section from a result bundle. The
%  output of the function handle must be either true or false, where true
%  indicates that the section is selected for extraction of results
%
%  Example function handle:
%
%  selection_function = @(section) section.video_properties.duration >= 10;
%  
%  This function selects only videos that have a duration of 10 seconds or
%  longer
    

%% Preliminaries

% ---
% Check if a selection function should be applied to eeach section

if nargin > 6 && ~isempty(varargin{1}) ...
        && strcmp(class(varargin{1}),'function_handle')
    % Apply a selection function to each section
    apply_selection = true;
    selection_function = varargin{1};
else
    % Do not apply selection function to sections
    apply_selection = false;
end

% ---
% Result bundle acquisition

if iscell(input_bundle)
    bundle_container = cell(1,numel(input_bundle));
    multiple_bundles = true;
    number_of_bundles = numel(input_bundle);
else
    bundle_container = cell(1,1);
    multiple_bundles = false;
    input_bundle = {input_bundle};
    number_of_bundles = 1;
end

for bb = 1:number_of_bundles
    
    this_bundle = input_bundle{bb};
    
    if isempty(this_bundle)
        % If an empty array is given as input_handle, choose an input bundle
        % from disk
        [load_file,load_path] = ...
            uigetfile('Select result bundle file from disk.');
        full_load_path = [load_path load_file];
        %load result bundle keyword table
        keyword_table = load(full_load_path,'keyword_table');
        keyword_table = keyword_table.keyword_table
        
    elseif ischar(this_bundle)
        % If a string is given as the input bundle, check if the string points
        % towards a directory, and if so let user pick a result bundle file
        % starting from that directory. If the string points towards a
        % .mat-file, try to use this file as a result bundle file to load
        is_mat_file = exist(this_bundle,'file');
        
        if is_mat_file==0 %if no file was found
            error('File %s does not exist.',this_bundle)
        elseif is_mat_file==7 %means that path is a directory
            % Pick file dialog with pre-set start directory
            [load_file,load_path] = ...
                uigetfile(this_bundle,'Select result bundle file from disk.');
            full_load_path = [load_path load_file];
            %load result bundle keyword table
            keyword_table = load(full_load_path,'keyword_table');
            keyword_table = keyword_table.keyword_table;
        elseif is_mat_file==2 && strcmp(this_bundle(end-3:end),'.mat')
            %means that path is a file whose path ends on '.mat'
            %load result bundle keyword table
            full_load_path = this_bundle;
            keyword_table = load(full_load_path,'keyword_table');
            keyword_table = keyword_table.keyword_table;
        end
    end
    
    
    % Number of sections in the bundle
    sections_in_bundle = numel(keyword_table)./2;
    
    % ---
    
    
    
    %% Apply keyword specifications to keyword_table and find conform sections
    
    % Start with all sections included
    query_inds = 1:sections_in_bundle;
    
    % ----
    % Apply sufficient_keywords to reject all sections that do not have any of
    % the keywords
    if ~isempty(sufficient_keywords)
        with_keyword_function = @(section) ...
            (sum(ismember(sufficient_keywords,section))>0);
        query_inds = intersect(query_inds,...
            find(cellfun(with_keyword_function,keyword_table(:,1))));
    end
    
    % ----
    % Apply necessary_keywords to reject all sections that do not have all
    % keywords
    if ~isempty(necessary_keywords)
        with_keyword_function = @(section) ...
            ~(sum(~ismember(necessary_keywords,section))>0);
        query_inds = intersect(query_inds,...
            find(cellfun(with_keyword_function,keyword_table(:,1))));
    end
    
    % ---
    % Apply not_keywords to reject all sections that contain any of the
    % keywords
    if ~isempty(not_keywords)
        % Rejection procedure
        with_keyword_function = @(section) ...
            sum(ismember(not_keywords,section))>0;
        include_not_keyword_inds = ...
            find(cellfun(with_keyword_function,keyword_table(:,1)));
        query_inds = setdiff(query_inds,include_not_keyword_inds);
    end
    
    % ---
    
    
    %% Execute the result query on each conform section
    
    number_returned = numel(query_inds);
    if number_returned == 0
        % If the query yielded no results from the keyword_table, assign
        % function value as false and stop further execution of function
        fprintf('No sections corresponding to keyword query were found.\n')
        fprintf('Query:\n')
        fprintf('OR keywords:\n')
        disp(sufficient_keywords)
        fprintf('AND keywords:\n')
        disp(necessary_keywords)
        fprintf('NOT keywords:\n')
        disp(not_keywords)
        bundle_container{bb} = false;
        break
    else
        % If the query yielded results from the keyword_table
        
        % ---
        % Check if Tracedropper masks should be considered, and if so, apply
        % the tracedrop_mask for each section
        
        if ~exist('non_dropped','var') || isempty(non_dropped) || ...
                ~non_dropped
            non_dropped = false;
        end
        
        % ---
        % Extract results
        
        % inclusion/rejection arrays based on existence of tracedrop mask
        % and conformity with a selection function
        mask_yesno = zeros(1,number_returned);
        selection_yesno = zeros(1,number_returned);
        
        % Find the sections that have a tracedrop_mask and that pass the
        % selection function
        for kk = 1:number_returned
            % ---
            % Check if section has tracedrop_mask as a field
            
            % Determine index of this section
            variable_to_load = keyword_table{query_inds(kk),2};
            % Load only this section from bundle
            this_section = ...
                load(full_load_path,variable_to_load);
            this_section = this_section.(variable_to_load);
            % Check if tracedrop_mask is actually a field
            mask_yesno(kk) = isfield(this_section,'tracedrop_mask');
            if apply_selection == true
                % Check if selection function returns true
                if selection_function(this_section)
                    selection_yesno(kk) = true;
                else
                    selection_yesno(kk) = false;
                end
            end
        end
        
        if non_dropped == false && apply_selection == false
            % Only Tracedropper application should be considered, but
            % a selection function is not applied
            
            % Keep only the indices to sections that do have a tracedrop
            % mask
            query_inds = query_inds(mask_yesno~=0);
            
        elseif non_dropped == true && apply_selection == true
            % TraceDropper application is not considered, but the selection
            % function is applied
            
            % Keep only the indices that have been selected by the
            % selection function
            query_inds = query_inds(selection_yesno);
            
        elseif non_dropped == false && apply_selection == true
            % TraceDropper is considered, and the selection function is
            % applied
            
            % Keep only the indices that have been selected by the
            % selection function and that have a tracedrop mask
            query_inds = query_inds(selection_yesno&mask_yesno~=0);
            
        end
        
        if isempty(query_inds) % How many sections are still in?
            % If none of the sections were selected, the function
            % returns with false. This means that no valid results are
            % contained in the query return
            bundle_container{bb} = false;
            break
        end
        
    end
    
    % Allocate cell array to contain the results for each section
    number_returned = numel(query_inds);
    this_extracted_results = cell(1,number_returned);
    
    
    % ----
    % Extract conform sections of the bundle using the query indices
    for kk = 1:number_returned
        % ---
        % Load the section from the saved bundle
        
        % Determine index of this section
        variable_to_load = keyword_table{query_inds(kk),2};
        % Load only this section from bundle
        this_section = ...
            load(full_load_path,variable_to_load);
        this_section = this_section.(variable_to_load);
        % Apply query function to get the requested results
        this_extracted_results{kk} = [feval(query_function,this_section)];
        if ~non_dropped % Apply tracedrop_mask to reject traces
            this_mask = this_section.tracedrop_mask;
            if numel(this_mask) == numel(this_extracted_results{kk})
                %If the extracted results are an array that the mask is
                %applicable to
                this_extracted_results{kk} = ...
                    this_extracted_results{kk}(this_section.tracedrop_mask);
            end
        end
    end
    bundle_container{bb} = this_extracted_results;
end

% Assign the extracted results from all sections into a cell array, one
% cell element per section
extracted_results = [bundle_container{:}];