function bundle_to_queried_bundle( ...
    input_bundle, target_bundle, ...
    sufficient_keywords,necessary_keywords,not_keywords)
%BUNDLE_TO_QUERIED_BUNDLE Save only parts of a result bundle according to a
%set of keywords
%
%   Detailed explanation goes here
%
%   Inputs:
%
%   input_path - path to the bundle which serves as the source, must be
%   as a string, if empty [] user is asked for input bundle, if only path
%   to a directory user will be asked to pick a file, but starting from
%   this directory
%
%   target_bundle - path to the bundle that will be saved to, must be
%   as a string, if empty [] user is asked for input bundle, if only path
%   to a directory user will be asked to pick a file, but starting from
%   this directory
%
%   sufficient_keywords
%   is a cell array of strings containing all keywords
%   that by themselves select a result section
%
%   necessary_keywords
%   is a cell array of strings containing all keywords
%   which a section must have to be extracted
%
%   not_keywords
%   is a cell array of strings that must not occur in the
%   keyword list



%% Preliminaries

% ---
% RESULT bundle acquisition

if isempty(input_bundle)
    % If an empty array is given as input_handle, choose an input bundle
    % from disk
    [load_file,load_path] = ...
        uigetfile('Select result bundle file from disk.');
    full_load_path = [load_path load_file];
    %load result bundle keyword table
    keyword_table = load(full_load_path,'keyword_table');
    keyword_table = keyword_table.keyword_table;
    
elseif ischar(input_bundle)
    % If a string is given as the input bundle, check if the string points
    % towards a directory, and if so let user pick a result bundle file
    % starting from that directory. If the string points towards a
    % .mat-file, try to use this file as a result bundle file to load
    is_mat_file = exist(input_bundle,'file');
    
    if is_mat_file==7 %means that path is a directory
        % Pick file dialog with pre-set start directory
        [load_file,load_path] = ...
            uigetfile(input_bundle,'Select result bundle file from disk.');
        full_load_path = [load_path load_file];
        %load result bundle keyword table
        keyword_table = load(full_load_path,'keyword_table');
        keyword_table = keyword_table.keyword_table;
    elseif is_mat_file==2 && strcmp(input_bundle(end-3:end),'.mat')
        %means that path is a file whose path ends on '.mat'
        %load result bundle keyword table
        full_load_path = input_bundle;
        keyword_table = load(full_load_path,'keyword_table');
        keyword_table = keyword_table.keyword_table;
    end
end

% Number of sections in the source bundle
sections_in_bundle = numel(keyword_table)./2;

% ---
% TARGET bundle acquisition

if isempty(target_bundle)
    % If an empty array is given as input_handle, choose an input bundle
    % from disk
    [target_file,target_path] = ...
        uigetfile('Where to save extracted bundle.');
    full_target_path = [target_path target_file];
    
elseif ischar(target_bundle)
    % If a string is given as the input bundle, check if the string points
    % towards a directory, and if so let user pick a result bundle file
    % starting from that directory. If the string points towards a
    % .mat-file, try to use this file as a result bundle file to load
    
    
    if strcmp(target_bundle(end-3:end),'.mat')
        %means that path is a file whose path ends on '.mat'
        %load result bundle keyword table
        full_target_path = target_bundle;
    else %means that path is a directory
        % Pick file dialog with pre-set start directory
        [target_file,target_path] = ...
            uigetfile('Where to save extracted bundle.');
        full_target_path = [target_path target_file];
    end
end

% ---

% Check if input bundle and target bundle are the same, if so, stop
% function
if strcmp(input_bundle,target_bundle)
    msgbox('Target and source bundle must be different!','Icon','Error')
    return
end
    




%% Apply keyword specifications to keyword_table and find conform sections

fprintf('Query based on keywords has been issued...\n')

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



%% Construct target bundle from source bundle


number_returned = numel(query_inds);


if number_returned == 0
    % If the query yielded no results from the keyword_table
    fprintf('No result sections were found for these keywords.\n')
else
    % If the query returned result sections, create target bundle containing
    % these result sections
    
    fprintf('%d result sections returned for query.\n',number_returned)
    
    % Store source keyword table in container variable
    source_keyword_table = keyword_table;
    
    % Allocate new keyword table for target bundle
    keyword_table = cell(number_returned,2);

    %Create new result bundle file, also clears out an existent bundle file
    %if it exists
    save(target_bundle,'keyword_table')
        
    % ----
    % Append returned result sections to target bundle
    for kk = 1:number_returned
        
        % --- Construct index under which the section will be saved in the
        % target bundle
        variable_to_store = ['section_' num2str(kk)];
        
        % ---
        % Load the section from the saved bundle
        variable_to_load = source_keyword_table{query_inds(kk),2};
        loaded_section= load(full_load_path,variable_to_load);
        eval([variable_to_store '= loaded_section.(variable_to_load);'])

        % ---
        % Append to target bundle
        
        % Store index in the keyword table
        keyword_table{kk,2} = variable_to_store;
        % Store keywords in the keyword table
        keyword_table{kk,1} = source_keyword_table{query_inds(kk),1};
        % Append variable
        save(target_bundle,variable_to_store,'-append');
        
        % Report progress of bundle creation
        fprintf('%d of %d sections added to bundle.\n',...
            kk,number_returned)
        
    end
    
end

% Save the keyword table to the file on disk that holds the bundled results
save(target_bundle,'keyword_table','-append')
fprintf('Creation of bundle %s finished.\n',full_target_path)