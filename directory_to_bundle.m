function ...
    directory_to_bundle(root_directory,target_file)
%DIRECTORY_TO_BUNDLE Bundles result files inside a directory tree into one
%result bundle
%   root_directory specifies a parent (or root) directory. All
%   Analysis_Results.mat files in this directory and its subdirectories are
%   found and treated as sections of data.
%
%   They are labelled with an increasing integer index
%   An index variable keyword_table is created that stores the keywords
%   associated with each integer index, and then saved to a file on disk
%   Each section is then saved as that integer index into the same file
%   saved on disk
%
%   target_directory specifies the path at which the results should be
%   saved


%% Find the paths of all Analysis_Results.mat files recursively

% If user did not supply inputs for root directory and/or target file,
% request them from user
if ~exist('root_directory','var')
    root_directory = uigetdir(pwd,'Choose root of directory tree.');
end
if ~exist('target_file','var')
    [target,target_path] = ...
        uiputfile(root_directory,'Choose file to save results to.');
    target_file = [target_path target];
end


% Initialization of cell array to contain the result paths
result_paths = {};

% Recursive exploration of directories

directory_batch = {root_directory};

fprintf('Recursive exploration of directory started...\n')

while ~isempty(directory_batch)
    
    % Draw last directory from directory cell
    this_directory = dir(directory_batch{end});
    this_directory_path = directory_batch{end};
    
    % Remove last directory from batch
    if numel(directory_batch)>1
        directory_batch = directory_batch(1:end-1);
    else
        directory_batch = {};
    end
    
    % Get subdirectories of this directory
    this_subdirectories = {this_directory([this_directory.isdir]).name};
    
    % Remove . and .. directories
    this_subdirectories = this_subdirectories( ...
        ~strcmp('.',this_subdirectories) & ...
        ~strcmp('..',this_subdirectories));
    
    if ~isempty(this_subdirectories)
        
        % Add path of this directory in front of the subdirectories
        this_path_addition = cell(1,numel(this_subdirectories));
        [this_path_addition{:}] = deal([this_directory_path filesep]);
        this_subdirectories = strcat(this_path_addition,this_subdirectories);
        
        % Add to the end of the directory batch
        directory_batch = {directory_batch{:},this_subdirectories{:}};
        
    end
    
    % Get paths of Analysis_Results.mat files from this directory
    this_result_paths = ...
        dir([this_directory_path filesep 'Analysis_Results.mat']);
    % If there is a result file(s), add it (them) to the list
    if ~isempty(this_result_paths)
        this_result_addition = cell(1,numel(this_result_paths));
        [this_result_addition{:}] = ...
            deal([this_directory_path filesep]);
        this_result_paths = strcat(this_result_addition, ...
            {this_result_paths.name});
        result_paths = {result_paths{:},this_result_paths{:}};
    end
    
end

% Create relative folder paths from root directory
root_path_length = length(root_directory);
rel_paths_function = @(full_path) full_path(root_path_length+2:end);
relative_paths = cellfun(rel_paths_function,result_paths,...
    'UniformOutput',false);

%Number of files that have been found
number_of_files = numel(relative_paths);

fprintf('%d "All_Results.mat" files found.\n',number_of_files)

%% Load result files, save into new result bundle file

% Allocate keyword table
keyword_table = cell(number_of_files,2);


%Create new result bundle file, also clears out an existent bundle file if
%it exists
save(target_file,'keyword_table')
fprintf('New bundle %s created, appending sections now...\n',target_file)


for ff = 1:number_of_files
    
    fprintf('Appending file %d of %d.\n',ff,number_of_files)
    
    index = ['section_' num2str(ff)];
    
    eval([index ' = struct;'])
    
    % Load the results for this file and store them into the variable
    % named according to the index ff
    this_result_struct = load(result_paths{ff});
    eval([index '.video_properties =' ...
        'this_result_struct.trace_results.video_properties;'])
    eval([index '.trace_results =' ...
        'this_result_struct.trace_results.trace_results;'])
    eval([index '.breakage_results =' ...
        'this_result_struct.breakage_results.breakage_results;'])
    if isfield(this_result_struct,'tracedrop_mask')
        eval([index '.tracedrop_mask =' ...
            'this_result_struct.tracedrop_mask;'])
    end
    
    % Add the variable to the file on disk that stores the bundled results
    save(target_file,index,'-append')
    
    % ---
    % Get and store keywords from file paths
    % Construct keyword string cell
    number_of_keywords = sum(relative_paths{ff}==filesep);
    keywords = textscan(relative_paths{ff},'%s', ...
        number_of_keywords, ...
        'Delimiter',filesep);
    % Store keywords in the keyword table
    keyword_table{ff,1} = keywords{:};
    % Store the variable name associated with these keywords
    keyword_table{ff,2} = index;
    
end

% Save the keyword table to the file on disk that holds the bundled results
save(target_file,'keyword_table','-append')
fprintf('Creation of bundle %s finished.\n',target_file)