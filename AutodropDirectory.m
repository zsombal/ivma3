% Which model should be used for tracedropping?
[file_name,file_path] = uigetfile(pwd,...
    'Select tracedrop training file...');
full_model_path = [file_path file_name];

% Which directory should be tracedropped?
parent_directory = uigetdir(file_path);

% Overwrite existent tracedrop_masks?
overwrite = 'Yes';

% ---
% Load existent model
training_expertise = load(full_model_path);

% ---
% Call function that automatically drops traces for all
% Analysis_Results.mat files in parent_directory and its subfolders

% Use optimal threshold in detection score that is saved with the model
% elapsed_time = ...
%     auto_tracedrop(training_expertise,parent_directory,overwrite);

% Use user-specified threshold for detection score
opt_thresh = 0.7;
elapsed_time = ...
    auto_tracedrop(training_expertise,parent_directory,overwrite,opt_thresh);

fprintf('Elapsed time: %f seconds.\n',elapsed_time);