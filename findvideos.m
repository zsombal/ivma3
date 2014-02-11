function [avi_file_paths,target_folder_paths] = ...
    findvideos( source_root_directory, target_root_directory )
%FINDVIDEOS Finds all the .avi videos in a parent folder and its
%subfolders
%   For an input directory source_root_directory (string) a cell containing the
%   paths to all .avi files in the parent folder specified by
%   source_root_directory is returned
%   target_folder_paths is a cell array containing the paths to the folders
%   that will contain the analysis results. They are assigned so as to
%   recreated the tree formed by all raw videos in source_root_directory

% Initialization of cell array to contain the .avi video paths
avi_file_paths = {};

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
    
    % Get paths of .avi files from this directory
    this_avi_file_paths = ...
        dir([this_directory_path filesep '*.avi']);
    % If there are .avi files, add them to the list
    if ~isempty(this_avi_file_paths)
        this_avi_addition = cell(1,numel(this_avi_file_paths));
        [this_avi_addition{:}] = ...
            deal([this_directory_path filesep]);
        this_avi_file_paths = strcat(this_avi_addition, ...
            {this_avi_file_paths.name});
        avi_file_paths = {avi_file_paths{:},this_avi_file_paths{:}};
    end
    
end

% Create relative folder paths from root directory
root_path_length = length(source_root_directory);
rel_paths_function = @(full_path) full_path(root_path_length+2:end);
avi_relative_paths = cellfun(rel_paths_function,avi_file_paths,...
    'UniformOutput',false);

% Remove all the paths to control videos that are produced by the
% EvaluateVideo function during analysis
control_string = '_Control.avi'; %String used to identify control videos
% Inverse control string for comparison operation
control_string = control_string(end:-1:1);
control_string_length = length(control_string);
% Function to detect trailing control video marker sequence
is_not_control_function = @(input_string) ...
    ~strncmp(control_string,input_string(end:-1:1),control_string_length);
% Detect indices of non-control videos
not_control_inds = cellfun(is_not_control_function,avi_relative_paths);
% Keep only those paths for non-control videos
avi_file_paths = avi_file_paths(not_control_inds);
target_folder_paths = avi_relative_paths(not_control_inds);



% Strip off .avi in the end of the relative path and put the target
% directory in front
stripped_paths_function = ...
    @(full_path) [target_root_directory filesep full_path(1:end-4)];
target_folder_paths = cellfun(stripped_paths_function,...
    target_folder_paths,...
    'UniformOutput',false);