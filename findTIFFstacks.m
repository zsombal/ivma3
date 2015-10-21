function [tiff_stack_paths,target_folder_paths] = ...
    findTIFFstacks( source_root_directory, target_root_directory )
%FINDTIFFSTACKS Finds all the .stk stacks in a parent folder and its
%subfolders
%   For an input directory source_root_directory (string) a cell containing the
%   paths to all .stk files in the parent folder specified by
%   source_root_directory is returned
%   target_folder_paths is a cell array containing the paths to the folders
%   that will contain the analysis results. They are assigned so as to
%   recreate the tree formed by all stacks in source_root_directory

% Initialization of cell array to contain the .stk stack paths
tiff_stack_paths = {};

% Recursive exploration of directories

directory_batch = {source_root_directory};

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
    
    % Get paths of .stk files from this directory
    this_tiff_stack_paths = ...
        dir([this_directory_path filesep '*.stk']);
    % If there are .stk files, add them to the list
    if ~isempty(this_tiff_stack_paths)
        this_stk_addition = cell(1,numel(this_tiff_stack_paths));
        [this_stk_addition{:}] = ...
            deal([this_directory_path filesep]);
        this_tiff_stack_paths = strcat(this_stk_addition, ...
            {this_tiff_stack_paths.name});
        tiff_stack_paths = {tiff_stack_paths{:},this_tiff_stack_paths{:}};
    end
    
end

% Create relative folder paths from root directory
root_path_length = length(source_root_directory);
rel_paths_function = @(full_path) full_path(root_path_length+2:end);
stk_relative_paths = cellfun(rel_paths_function,tiff_stack_paths,...
    'UniformOutput',false);

% Create paths to the target folders
target_folder_paths = stk_relative_paths;


% Strip off .stk in the end of the relative path and put the target
% directory in front
stripped_paths_function = ...
    @(full_path) [target_root_directory filesep full_path(1:end-4)];
target_folder_paths = cellfun(stripped_paths_function,...
    target_folder_paths,...
    'UniformOutput',false);