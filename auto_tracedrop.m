function [ elapsed_time ] = auto_tracedrop( ...
    training_expertise,parent_directory,overwrite,opt_thresh)
%AUTO_TRACEDROP Uses a pretrained bagged decision tree model to drop traces
%
% This script will address all results (that have not yet been bundled) by
% adding tracedrop mask to them. A tracedrop mask contains the information
% if a specific filament's trace looks "proper" or if it is indicative of
% irregular movement of the filament during observation of the filament.
%
% Inputs:
%
% training_expertise (structure array)
% A training expertise file that has been created before by training on
% result bundles created from data scored by a human operator. The actual
% variable is passed, not the path string.
%
% parent_directory (path string)
% The script will address all result files contained in the parent
% directory and all its subdirectories (recursive search).
%
% overwrite ( string, "Yes" or "No" )
% Determines if existent tracedrop masks will be overwritten (Yes) or kept
% as they are (No).

% Preallocations
elapsed_time=-1;

% ---
% Retrieve data for classification from trained model
relevant_pixels = training_expertise.relevant_pixels;
if ~exist('opt_thresh','var')
    opt_thresh = training_expertise.opt_thresh;
end
decision_tree_model = training_expertise.decision_tree_model;
height = training_expertise.height;
width = training_expertise.width;





% ---
% How should existent tracedropper_masks be treated?
if strcmp(overwrite,'ask')
    % Ask user via dialogbox
    overwrite = questdlg('Should existent TraceDropper masks be replaced?',...
        'Yes','No');
end

tic

% ---
% Detect all Analysis_Results.mat files

% Initialization of cell array to contain the .avi video paths
file_paths = {};

% Recursive exploration of directories

directory_batch = {parent_directory};

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
    this_file_path = ...
        dir([this_directory_path filesep 'Analysis_Results.mat']);
    % If there is an Analysis_Results.mat files, add it to the list
    if ~isempty(this_file_path)
        file_paths = {file_paths{:},...
            [this_directory_path filesep 'Analysis_Results.mat']};
    end
    
end

number_of_files = numel(file_paths);

if number_of_files == 0
    % If no files were found, stop function execution
    fprintf('No Analysis_Results.mat file found.')
    msgbox('No Analysis_Results.mat file found.',...
        'No files found.','error')
    return
end


% ---
% Check for preexistent tracedrop_masks and videos of wrong dimension




fprintf('Checking for existence of tracedrop_masks and proper video dimensions...\n')

% Only result sections without tracedrop_mask are assigned a new
% tracedrop_mask
with_tracedrop_mask = false(1,number_of_files);
right_measure_inds = false(1,number_of_files);
for kk = 1:number_of_files
    results = ...
        (load(file_paths{kk}));
    with_tracedrop_mask(kk) = isfield(results,'tracedrop_mask');
    right_measure_inds(kk) = ...
        (results.trace_results.video_properties.width == width) && ...
        (results.trace_results.video_properties.height == height);
end

if ~strcmp(overwrite,'Yes')
    % Only keep those without tracedrop_mask and with proper dimensions
    file_paths = file_paths(~with_tracedrop_mask & right_measure_inds);
else
    % Only keep those with proper dimensions
    file_paths = file_paths(right_measure_inds);
end
old_number_of_files = number_of_files;
number_of_files = numel(file_paths);

if number_of_files == 0
    % If no files are left, stop function execution
    if strcmp(overwrite,'No')
        fprintf('\n%d of %d Analysis_Results.mat files already have tracedrop_masks\n',...
            sum(tracedrop_masks),old_number_of_files)
        fprintf('You chose not to overwrite any existent tracedrop_masks.\n')
    end
    fprintf('\n%d of %d Analysis_Results.mat files have wrong video dimension\n',...
        sum(right_measure_inds),old_number_of_files)
    fprintf('Automated dropping aborted.\n')
    msgbox(...
        ['All Analysis_Results.mat are either already tracedropped',...
        'or are of wrong video size. Check command line.'],...
        'No valid files.',...
        'error')
    return
end



fprintf('done.\n')

clear results



% ---
% Classify traces and store according tracedrop_mask

for kk = 1:number_of_files
    
    % ---
    % Extract/calculate features for this result file
    fprintf('Extracting trace features for file %d of %d.\n',...
        kk,number_of_files)
    
    if ~tracedrop_masks(kk) || strcmp(overwrite,'Yes')
        
        % If (1) no tracedrop mask exists OR (2) tracedrop mask should be
        % overwritten, go ahead and determine the tracedrop mask from the
        % training expertise, and place it with the results for this
        % specific video.
    
        % Load results from this file
        this_results = load(file_paths{kk});
        trace_results = this_results.trace_results.trace_results;

        % Extract tracedrop masks and trace image data and bounding boxes from
        % results bundle

        trace_images = {trace_results.trace_image};
        trace_bounding_boxes = {trace_results.trace_bounding_box};
        trace_lengths = [trace_results.trace_length];
        trace_widths = [trace_results.trace_width];
        trace_solidities = [trace_results.trace_solidity];

        % ---
        %Specifications for corner detection
        filter_coefficients = fspecial('gaussian',[21 1],2.5);
        maximal_corners = 200;

        % Determine number of featuresout-of-bundle data from
        number_of_features = 2.*maximal_corners + 6;

        number_of_traces = numel(trace_images);

        % Extract and/or store trace features
        feature_container = zeros(number_of_traces,number_of_features);
        % Create trace images that are centered and have same
        centered_images = false(number_of_traces,height.*width);

        parfor tt = 1:number_of_traces;


            % Initialize all 0 image
            this_image = false(height,width);

            % horizontal positioning
            xx_center = floor(trace_bounding_boxes{tt}(1) + ...
                trace_bounding_boxes{tt}(3)./2);
            xx_min = ceil(trace_bounding_boxes{tt}(1))-xx_center+round(width./2);
            xx_max = xx_min + trace_bounding_boxes{tt}(3)-1;
            % vertical positioning
            yy_center = floor(trace_bounding_boxes{tt}(2) + ...
                trace_bounding_boxes{tt}(4)./2);
            yy_min = ceil(trace_bounding_boxes{tt}(2))-yy_center+round(height./2);
            yy_max = yy_min + trace_bounding_boxes{tt}(4)-1;
            this_image(yy_min:yy_max,xx_min:xx_max) = ...
                trace_images{tt};

            % Get region properties of the trace in this image
            properties = regionprops(this_image,'Area','Orientation','Perimeter');

            % Rotate this image to have the trace main axis in horizontal direction
            this_image = imrotate(...
                this_image,-properties.Orientation,'bilinear','crop')

            %Store pixel data for this trace
            centered_images(tt,:) = this_image(:); %Store as row in centered_images

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

            feature_container(tt,:) = this_trace_features;

            % Plot the image of this trace with detected corners
            %     imagesc(this_image)
            %     hold on
            %     plot(corner_positions(:,1),corner_positions(:,2),'ro')
            %     hold off
            %     pause(1)



            fprintf('%d of %d traces done, file %d of %d.\n',...
            tt,number_of_traces,kk,number_of_files)

        end

        feature_container(:,2.*maximal_corners+2) = trace_lengths;
        feature_container(:,2.*maximal_corners+3) = trace_widths;
        feature_container(:,2.*maximal_corners+4) = trace_solidities;

        %Remove all imaginary parts from feature values
        feature_container = real(feature_container);

        % Clear data out of memory
        clear trace_images trace_bounding_boxes heights widths

        cropped_pixel_data = double(centered_images(:,relevant_pixels));

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
            predict(decision_tree_model,feature_container);
        predicted_class = logical(cellfun(@(xx) str2num(xx),predicted_class));
        fprintf('done.\n')

        % Determine tracedrop_mask according to classification
        tracedrop_mask = prediction_score(:,2)>=opt_thresh;

        % Append tracedrop_mask
        save(file_paths{kk},'tracedrop_mask','-append')
        
    end
    
end

elapsed_time = toc;