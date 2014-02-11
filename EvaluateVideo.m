function [ target_folder_list, error_log, varargout ] = EvaluateVideo(...
    videos,parameters,options)
%EVALUATEVIDEO Analyzes one or a batch of in vitro motility assay videos
%
%   INPUTS:
%
%   Videos: N-by-2 cell array, each of the N rows containing a video to be
%   analyzed and the target folder to save the analysis results to
%
%   Parameters: 1-by-1 or N-by-1 structure array, contains the analysis
%   parameters that are applied for all videos (1-by-1 case) or for each
%   individual video (N-by-1 case)
%
%   Options: 1-by-1 or N-by-1 structure array, contains options for saving
%   of intermediate control data during analysis
%
%   OUTPUTS:
%
%   target_folder_list: contains the paths to the folders into which
%   analysis results were saved successfully, N_s-by-1 cell array, where
%   N_s is the number of videos out of the whole batch, which could be
%   analyzed successfully
%
%   error_log: returns false if no errors have been logged, otherwise
%   N-by-1 cell array with strings containing error logs for each video in
%   the batch
%
%   Optional outputs - can consume a lot of memory when many videos are
%   included in the batch.
%
%   video_results_batch: N-by-1 cell array, each cell element contains the
%   filament motion analysis results for one video in teh form of a
%   structure array
%
%   breakage_results_batch: N-by-1 cell array, each cell element contains
%   the breakage analysis results for one video in the form of a cell array
%
%   elapsed_time: scaler double indicating the time taken for the complete
%   analysis of all videos, in seconds


%% Start clock for runtime estimation
tic;

%% Check if input arrays are right type and dimension

%Check if the videos batch conforms with the requirements
video_batch_size = size(videos);
if numel(size(video_batch_size))~=2 || video_batch_size(2) ~= 2 || ~(iscell(videos))
    error('ivma3:EvaluateVideo:ImproperVideoBatch',...
        'The video batch must be 2-by-N cell array for N videos.')
end

%Read parameters from parameter structure array
parameter_size = size(parameters);
if ~isstruct(parameters) || parameter_size(2)~=1
    error('ivma3:EvaluateVideo:ImproperParameterStructureArray',...
        'The parameters must be passed as a 1-by-N structure array.')
end

%Read options from option structure array
option_size = size(options);
if ~isstruct(options) || option_size(2)~=1
    error('ivma3:EvaluateVideo:ImproperOptionStructureArray',...
        'The parameters must be passed as a 1-by-N structure array.')
end

%Check if parameter and option structure array apply for all videos, or if
%specific values for each video are present. Otherwise return an error.
if parameter_size(1)==1
    general_parameters = true; %Apply same parameters to all videos
elseif parameter_size(1)==video_batch_size(1)
    general_parameters = false; %Apply different parameters for each video
else
    error('ivma3:EvaluateVideo:ImproperParameterStructureArraySize',...
        'Parameter array must be 1-by-1 or N-by-1, N=number of videos.')
end
if option_size(1)==1
    general_options = true; %Apply same options to all videos
elseif option_size(1)==video_batch_size(1)
    general_options = false; %Apply different options for each video
else
    error('ivma3:EvaluateVideo:ImproperOptionStructureArraySize',...
        'Option array must be 1-by-1 or N-by-1, N=number of videos.')
end

%Number of videos in the batch
videos_in_batch = video_batch_size(1);

%% Determine which outputs are requested by the function call

trace_result_batch_requested = false;
breakage_result_batch_requested = false;
elapsed_time_requested = false;
number_of_outputs = nargout;
if number_of_outputs > 2
    trace_result_batch_requested = true;
end
if number_of_outputs > 3
    breakage_result_batch_requested = true;
end
if number_of_outputs > 4
    elapsed_time_requested = true;
end


%% Array allocation

%To contain error message strings for each video in the batch
error_log = cell(1,video_batch_size(1));
% Error flag for each video
no_error_flag = true(1,video_batch_size(1));

%To contain the results for each video in the batch
if trace_result_batch_requested
    trace_results_batch = cell(1,video_batch_size(1));
end
if breakage_result_batch_requested
    breakage_results_batch = cell(1,video_batch_size(1));
end

parfor vv = 1:videos_in_batch
    
    try
        %Encloses the following operations to catch errors instead of
        %interrupting the processing of the whole batch
        
        %% Initialize switches etc. for this loop iteration
        
        % Switch to tell if immotile frame removal was successful
        immotile_removal_success = false;
        % Thrshold velocity to separate immotile frames from motile frames
        imm_threshold_velocity = 0;
        % Proportion of frame-to-frame velocities that are above the motility
        % threshold
        motile_fraction = 0;
        % Video that is analyzed in this loop
        current_video = [];
        % v_max in case immotile frame removal is activated
        frame_to_frame_vmax = [];
        
        
        %% Read analysis parameters from input variables
        
        % Find out if general parameters (applying to all videos in the same way)
        % or specific parameter sets (specific for each video) are passed
        % into the function
        if general_parameters==false
            param_ind = vv;
        else
            param_ind = 1;
        end
        
        %Micrometers per pixel in the video
        experimental_scaling_factor = ...
            parameters(param_ind).scaling_factor;
        
        %Raw frames merged into one compressed frame
        frames_to_merge = ...
            parameters(param_ind).frames_to_merge;
        
        %Black-white threshold value, between 0 and 1 (black-to-white
        %grey-scales, respectively)
        BW_threshold = ...
            parameters(param_ind).BW_threshold;
        
        %Minimal length below which detected objects are excluded from the
        %analysis (length in micrometers)
        min_length = ...
            parameters(param_ind).min_length;
        
        %Maximal relative length change up to which filaments are tracked as
        %the same filament (length in micrometers)
        rel_length_change = ...
            parameters(param_ind).rel_length_change;
        
        %Maximal absolute length change up to which filaments are tracked as
        %the same filament (length in micrometers)
        abs_length_change = ...
            parameters(param_ind).abs_length_change;
        
        %Minimum time (seconds) a filament has to be present to be included
        %for analysis
        min_presence_time = ...
            parameters(param_ind).min_presence_time;
        
        %Minimal length of a trace to be included for analysis
        %(length in micrometers)
        min_trace_length = ...
            parameters(param_ind).min_trace_length;
        
        %Removal of immotile frames ON (true) or OFF (false)
        immotile_frame_removal = ...
            parameters(param_ind).immotile_frame_removal;
        
        % Find out if general options (applying to all videos)
        % or specific options (for each video) sets are passed
        % into the function
        if general_options==false
            opt_ind = vv;
        else
            opt_ind = 1;
        end
        
        %Should a control video of the preprocessing be created in the target
        %folder?
        preprocessing_control_video = options(opt_ind).preprocessing_control_video;
        
        %Should a control video of filament tracking be created in the target
        %folder?
        tracking_control_video = options(opt_ind).tracking_control_video;
        
        %Should the breakage and trace results be stored in a plain text file?
        write_plain_text = options(opt_ind).write_plain_text;
        
        
        %% Opening of video file as VideoReader object
        
        do_analysis = true;
        try %Try if video can be read
            current_video = VideoReader(videos{vv,1});
        catch ME %otherwise store error message and prevsnt further analysis
            error_log{vv} = sprintf('Error accessing video file:\n%s',videos{vv,1});
            do_analysis = false;
        end
        
        if do_analysis
            
            %% Creation of target directory to save to
            this_target_folder = videos{vv,2};
            [~,~,~] = mkdir(this_target_folder);
            
            
            %% Store analysis parameters in plain text file
            fid = fopen([this_target_folder filesep 'Parameters.txt'],'w');
            this_video_parameter_fields = fieldnames(parameters(param_ind));
            fprintf(fid,'Input video:\n%s\n\n',videos{vv,1});
            fprintf(fid,'\nAnalysis parameters:\n');
            for field = this_video_parameter_fields.'
                fprintf(fid,'%s = %f\n',...
                    char(field),...
                    parameters(param_ind).(char(field)));
            end
            fclose(fid);
            
            
            %% Acquisition of frames from VideoReader object
            
            % --------------
            % Video specs
            %Get the number of frames, frame rate, and total elapsed time
            number_of_frames = current_video.NumberOfFrames;
            frame_rate = current_video.FrameRate;
            duration = current_video.Duration;
            %Video size
            [width,height] = deal(current_video.Width,...
                current_video.Height);
            
            
            %% Frame merging, brightness balancing, rescaling
            
            number_of_merged_frames = floor(number_of_frames./frames_to_merge)-1;
            % Stop analyzing this video and write error report if only one
            % merged frame is left.
            if number_of_merged_frames < 2
                no_error_flag(vv) = false;
                error_log{vv} = sprintf(...
                    ['Not enough frames for analysis, only %d frames after merging, in\n%s\n' ...
                    'Decrease frames_to_merge number or take longer video.'], ...
                    number_of_merged_frames,[this_target_folder]);
                fid = fopen([this_target_folder filesep 'error_report.txt'],'w');
                fprintf(fid,'%s\n',error_log{vv});
                fclose(fid);
                continue
            end
            
            merged_frames = zeros(height,width,number_of_merged_frames);
            merged_frame_brightness_medians = zeros(1,number_of_merged_frames);
            
            % ---
            %Merge frames and get brightness medians
            for mm = 1:number_of_merged_frames
                
                % -------------
                % Read video frames into matrix
                % Matrix containing all frame data for frames merged into this
                % merged frame, dimensions are:
                % image-height,image-width,RGB-channel,frame
                % Take mean of the three color channels to get gray scale frame
                % data, also scale from [0,255] to [0,1] range
                current_gray_frames = ...
                    mean( ...
                    read(current_video, ...
                    [(mm.*frames_to_merge),((mm+1).*frames_to_merge-1)]),3)./255
                
                %Take the mean brightness of the frames to be merged and assign
                %it to the merged frame
                merged_frames(:,:,mm) = ...
                    mean(squeeze(current_gray_frames(:,:,1,:)),3);
                brightnesses_in_this_merged_frame = merged_frames(:,:,mm);
                merged_frame_brightness_medians(mm) = ...
                    median(brightnesses_in_this_merged_frame(:));
            end
            
            % Remove current_video from memory
            current_video = [];
            
            % Remove gray scale array from memory
            current_gray_frames = [];
            
            mean_median_brightness = mean(merged_frame_brightness_medians);
            
            % ----
            %Shift to have the same brightness median, compensates for
            %photo-bleaching
            for mm = 1:number_of_merged_frames
                merged_frames(:,:,mm) = merged_frames(:,:,mm) ...
                    - (merged_frame_brightness_medians(mm) ...
                    -mean_median_brightness);
            end
            
            % ----
            %Rescaling so that the brightness range of the video fits exactly
            %into the available brightness range of [0,1]
            min_brightness = min(merged_frames(:));%Minimum of all brightnesses
            max_brightness = max(merged_frames(:));%Maximum of all brightnesses
            merged_frames(:) = (merged_frames(:)-min_brightness) ...
                ./ (max_brightness-min_brightness);%Rescaling
            
            % ----
            % Create a control video of the merging, rescaling
            % operations, useful for debugging
            if preprocessing_control_video
                
                try %write the preprocessing control video
                    
                    % Control output of rescaled video)
                    preprocessed_frames = reshape(merged_frames,...
                        height,width,1,number_of_merged_frames);
                    
                    % ----
                    % Open VideoWrite object
                    preprocessing_video = VideoWriter(...
                        [this_target_folder filesep 'PreProcessing_Control.avi'],...
                        'Uncompressed AVI');
                    preprocessing_video.FrameRate = frame_rate./frames_to_merge;
                    
                    % ----
                    %Open, write into, and close file associated with video object
                    open(preprocessing_video);
                    %preprocesses_frames has to be scaled to [0,255] range and then
                    %converted to unsigned integer to fit the uncompressed .avi
                    %format for output
                    writeVideo(preprocessing_video,uint8(preprocessed_frames.*255));
                    close(preprocessing_video)
                    % Remove preprocessing_video and preprocessed_frames
                    % from memory
                    preprocessing_video = [];
                    preprocessed_frames = [];
                    
                catch ME %make an entry into the error log
                    error_log{vv} = sprintf(...
                        'Error writing control video file :\n%s', ...
                        [this_target_folder filesep 'PreProcessing_Control.avi']);
                    no_error_flag(vv) = false;
                    fid = fopen([this_target_folder filesep 'error_report.txt'],'w');
                    fprintf(fid,'%s\n',error_log{vv});
                    fclose(fid);
                    continue %skip to next video in the batch
                end
                
            end
            
            
            %% Object detection, length and width calculation, length filtering
            
            % ----
            % Thresholding into frames containing binary images
            % false corresponds to no object, true corresponds to object
            binary_frames = merged_frames>BW_threshold;
            
            % ----
            %Inline functions to calculate object length and width from object
            %area and perimeter based on a transformation into a rectangle with
            %same area and perimeter
            rect_length = @(area,perim) ...
                experimental_scaling_factor.*(perim./4 + sqrt(perim.^2./16-area));
            rect_width = @(area,perim) ...
                experimental_scaling_factor.*(perim./4 - sqrt(perim.^2./16-area));
            
            % ----
            %Detection of contiguous true (white) regions aka connected
            %components and their properties
            
            %Numeric array to count the number of filaments in each frame
            filaments_in_frame = zeros(1,number_of_merged_frames,'uint32');
            
            %Cell array that will contain the properties for each of the merged
            %frames
            frame_conncomp_props = cell(1,number_of_merged_frames);
            
            %Structure array that will contain the results needed for
            %breakage analysis
            breakage_results = struct;
            breakage_results(number_of_merged_frames).times = [];
            %Store merged frame times for breakage analysis, in seconds
            [breakage_results.times] = deal((1:number_of_merged_frames-0.5) ...
                .*frames_to_merge./frame_rate);
            
            for ff = 1:number_of_merged_frames
                conn_comps = bwconncomp(binary_frames(:,:,ff),4);
                frame_conncomp_props{ff} = regionprops(conn_comps,...
                    merged_frames(:,:,ff),...
                    {'Area','Perimeter','WeightedCentroid', 'BoundingBox', ...
                    'PixelIdxList','PixelList','Image','Solidity'});
                
                try
                    %The 'try' is used here so that an error can be raised
                    %wherever all filaments are lost due to rejection of
                    %filaments
                    
                    if numel(frame_conncomp_props)==0
                        error('Out of object error')
                    end
                    
                    
                    %Store raw pixel counts of objects for breakage analysis
                    breakage_results(ff).areas = ...
                        [frame_conncomp_props{ff}.Area];
                    
                    % Keep only objects that are greater than one pixel
                    frame_conncomp_props{ff} = ...
                        frame_conncomp_props{ff}( ...
                        [frame_conncomp_props{ff}.Area]>1);
                    
                    
                    if numel(frame_conncomp_props)==0
                        error('Out of object error')
                    end
                    
                    %Calculate lengths and width of objects
                    these_areas = [frame_conncomp_props{ff}.Area];
                    these_perims = [frame_conncomp_props{ff}.Perimeter];
                    length_cell = num2cell(...
                        arrayfun(rect_length,these_areas,these_perims));
                    width_cell = num2cell(...
                        arrayfun(rect_width,these_areas,these_perims));
                    [frame_conncomp_props{ff}.Length] = ...
                        length_cell{:};
                    [frame_conncomp_props{ff}.Width] = ...
                        width_cell{:};
                    
                    %Store lengths and widths for breakage results
                    breakage_results(ff).lengths = [frame_conncomp_props{ff}.Length];
                    breakage_results(ff).widths = [frame_conncomp_props{ff}.Width];
                    
                    % Keep only those filaments whose length is equal or above a
                    % minimum filament length
                    [frame_conncomp_props{ff}] = ...
                        frame_conncomp_props{ff}( ...
                        [frame_conncomp_props{ff}.Length]>=min_length);
                    
                    
                    if numel(frame_conncomp_props)==0
                        error('Out of object error')
                    end
                    
                    % ----
                    % Sort out objects touching the edge of the video frame
                    
                    % Conversion from linear indices to two-dimensional indices
                    to_xy_conversion = @(lin_inds) ind2sub([height,width],lin_inds);
                    % Conversion for all filaments in this frame
                    [xx_inds,yy_inds] = ...
                        cellfun(to_xy_conversion,...
                        {frame_conncomp_props{ff}.PixelIdxList},...
                        'UniformOutput',false);
                    % Function returns true if one or more pixels are on the edge,
                    % and false if none are on the edge
                    on_the_edge = @(yy,xx) ...
                        sum(ismember(1,xx) | ismember(width,xx) ...
                        | ismember(2,xx) | ismember(width-1,xx) ...
                        | ismember(3,xx) | ismember(width-2,xx) ...
                        | ismember(4,xx) | ismember(width-3,xx) ...
                        | ismember(5,xx) | ismember(width-4,xx) ...
                        | ismember(6,xx) | ismember(width-5,xx) ...
                        | ismember(7,xx) | ismember(width-6,xx) ...
                        | ismember(8,xx) | ismember(width-7,xx) ...
                        | ismember(9,xx) | ismember(width-8,xx) ...
                        | ismember(10,xx) | ismember(width-9,xx) ...
                        | ismember(1,yy) | ismember(height,yy) ...
                        | ismember(2,yy) | ismember(height-1,yy) ...
                        | ismember(3,yy) | ismember(height-2,yy) ...
                        | ismember(4,yy) | ismember(height-3,yy) ...
                        | ismember(5,yy) | ismember(height-4,yy) ...
                        | ismember(6,yy) | ismember(height-5,yy) ...
                        | ismember(7,yy) | ismember(height-6,yy) ...
                        | ismember(8,yy) | ismember(height-7,yy) ...
                        | ismember(9,yy) | ismember(height-8,yy) ...
                        | ismember(10,yy) | ismember(height-9,yy) ...
                        );
                    
                    % Apply to all filaments in this frame
                    pixels_on_edge = cellfun(on_the_edge,xx_inds,yy_inds);
                    % Keep only those filaments that have no pixels touching the
                    % edge
                    none_on_the_edge_inds = find(pixels_on_edge==0);
                    frame_conncomp_props{ff} = ...
                        frame_conncomp_props{ff}(none_on_the_edge_inds);
                    
                    
                    if numel(frame_conncomp_props)==0
                        error('Out of object error')
                    end
                    
                    filaments_in_frame(ff) = ...
                        numel(none_on_the_edge_inds);
                    
                catch ME
                    
                    %This is executed if no filament was left after the above
                    %rejection procedures
                    
                    filaments_in_frame(ff) = 0;
                    
                end
                
            end
            
            % Clear merged frame array from memory
            merged_frames = [];
            
            %% Filament tracking with optional control output of marked video
            
            % ----
            % Cell array to contain the index forwarding vectors
            reverse_pointer_vectors = cell(1,number_of_merged_frames);
            
            %Counter to contain the number of all newly appearing filaments
            number_of_nonconnected_filaments = 0;
            
            % ----
            % Determine the dissimilarity matrix between each two consecutive
            % frames, infer a reverse pointer vector connecting filaments in
            % consecutive frames
            for tt = 1:number_of_merged_frames-1
                
                if filaments_in_frame(tt+1) == 0
                    %In case there are no filaments in the next frame, no
                    %reverse pointing vector is needed:
                    reverse_pointer_vectors{tt} = [];
                elseif filaments_in_frame(tt) == 0
                    % In case the last frame does not hold any filaments, all
                    % filaments in the next frame need a new index. Therefore
                    % all elements in the reverse pointer vector have to be -1
                    reverse_pointer_vectors{tt} = ...
                        -ones(filaments_in_frame(tt+1),1);
                else
                    %When both the last and the next frame contain filaments,
                    %the forwarding vector has to be assigned according to the
                    %dissimilarity matrix
                    
                    %Meshgrids of areas and centroid positions
                    [source_area_mesh,target_area_mesh] = ...
                        meshgrid([frame_conncomp_props{tt}.Area], ...
                        [frame_conncomp_props{tt+1}.Area]);
                    
                    source_frame_centroids = ...
                        [frame_conncomp_props{tt}.WeightedCentroid];
                    %Must be reshaped from 2N-by-1 array into N-by-2 array
                    source_frame_centroids = ...
                        [source_frame_centroids(1:2:end-1);...
                        source_frame_centroids(2:2:end)].';
                    
                    target_frame_centroids = ...
                        [frame_conncomp_props{tt+1}.WeightedCentroid];
                    %Must be reshaped from 2N-by-1 array into N-by-2 array
                    target_frame_centroids = ...
                        [target_frame_centroids(1:2:end-1);...
                        target_frame_centroids(2:2:end)].';
                    
                    
                    [source_centroidxx_mesh,target_centroidxx_mesh] = ...
                        meshgrid(source_frame_centroids(:,1),...
                        target_frame_centroids(:,1));
                    [source_centroidyy_mesh,target_centroidyy_mesh] = ...
                        meshgrid(source_frame_centroids(:,2),...
                        target_frame_centroids(:,2));
                    
                    % Matrix that assesses the dissimilarity between each
                    % combination of filaments in this and filaments in the next
                    % frame
                    dissimilarity_matrix = ...
                        (source_area_mesh-target_area_mesh).^2 ...
                        + (source_centroidxx_mesh-target_centroidxx_mesh).^2 ...
                        + (source_centroidyy_mesh-target_centroidyy_mesh).^2;
                    
                    % Matrix that contains -1 except for the column-wise minima, which
                    % hold the dissimilarity-value. This determines which filament
                    % in the next frame the filaments in this frame forward to
                    column_minima_matrix = -ones(size(dissimilarity_matrix));
                    [column_min_values,column_min_indices] = ...
                        min(dissimilarity_matrix,[],1);
                    column_minima_lin_inds = sub2ind(size(column_minima_matrix),...
                        column_min_indices,double([1:filaments_in_frame(tt)]));
                    column_minima_matrix(column_minima_lin_inds) = ...
                        column_min_values;
                    
                    % Vector that contains the indices of the minimum error
                    % filament from the last frame for the filaments in the next
                    % frame. Where no filament from the last frame points to the
                    % next filament, -1 is in the vector. The vector index
                    % represents filaments in the next frame, the content is the
                    % filament from the last frame that the index should be
                    % inherited from
                    
                    %Find minimal non-negative values and their indices. The
                    %inverse is taken to turn the minimal non-negative value into
                    %maxima, which are much easier to find.
                    [row_max_values,reverse_pointer_vectors{tt}] = ...
                        max(1./column_minima_matrix,[],2);
                    reverse_pointer_vectors{tt}(row_max_values==-1) = -1;
                    
                    % Apply the length change absolute and relative criterion to
                    % check if two connected filaments should actually be
                    % disconnected. If they are to different to be understood as
                    % the same filament in two consecutive frames, assign the
                    % forward vector as -1
                    
                    %Get the lengths
                    valid_reverse_pointer_indices = ...
                        find(reverse_pointer_vectors{tt}~=-1);
                    last_frame_lengths = ...
                        [frame_conncomp_props{tt}( ...
                        reverse_pointer_vectors{tt}(valid_reverse_pointer_indices) ...
                        ).Length];
                    next_frame_lengths = ...
                        [frame_conncomp_props{tt+1}( ...
                        valid_reverse_pointer_indices ...
                        ).Length];
                    
                    abs_length_diffs = abs(next_frame_lengths-last_frame_lengths);
                    rel_length_diffs = abs_length_diffs./ ...
                        (last_frame_lengths+next_frame_lengths).*2;
                    
                    % Apply the overlap criterion - two filaments that are about
                    % to be connected between the last and the next frame that do
                    % NOT have an overlap are disconnected. This means that their
                    % position in the forwarding vector is assigned -1
                    
                    % Get cell arrays containing the linear indices of filaments in
                    % the last and in the next frame
                    
                    last_frame_indices = {frame_conncomp_props{tt}( ...
                        reverse_pointer_vectors{tt}(valid_reverse_pointer_indices) ...
                        ).PixelIdxList};
                    next_frame_indices = {frame_conncomp_props{tt+1}( ...
                        valid_reverse_pointer_indices ...
                        ).PixelIdxList};
                    overlaps = false(1,numel(valid_reverse_pointer_indices));
                    for oo = 1:numel(overlaps)
                        common_pixels = intersect(last_frame_indices{oo}, ...
                            next_frame_indices{oo});
                        overlaps(oo) = numel(common_pixels)>0;
                    end
                    
                    %Assign the reverse pointer elements, that connect filaments
                    %which are too different, with the value -1 to disconnect the
                    %filaments
                    reverse_pointer_vectors{tt}( valid_reverse_pointer_indices( ...
                        abs_length_diffs>abs_length_change ...
                        | rel_length_diffs>rel_length_change ...
                        | 1./rel_length_diffs<rel_length_change ...
                        | ~overlaps)) ...
                        = -1;
                    
                end
                
                this_frame_nonconnected_filaments = ...
                    find(reverse_pointer_vectors{tt}==-1);
                number_of_nonconnected_filaments = ...
                    number_of_nonconnected_filaments ...
                    + numel(this_frame_nonconnected_filaments);
                
                
                
            end
            
            %% Connect tracked filaments and merge into traces
            
            % ----
            %Cell array to contain the tracked filament indices in each of the
            %frames. Each cell represents a frame of the video, and inside
            %there, the filaments are assigned an index. If the same index
            %appears in two different frames, it means that the objects in these
            %two frames, addressed by that index, are connected as the same
            %object, only observed in different frame.
            
            tracked_indices_array = cell(1,number_of_merged_frames);
            % To keep track of assigned indices
            current_index_counter = filaments_in_frame(1);
            tracked_indices_array{1} = [1:filaments_in_frame(1)];
            
            for tt = 1:number_of_merged_frames-1
                
                % ---
                %First, assign the forwarded indices
                
                % Find valid forwarding indices
                valid_reverse_vector_inds = ...
                    find(reverse_pointer_vectors{tt}~=-1);
                % Find which filaments in the last should be forwarded from
                valid_reverse_vector_pointers = ...
                    reverse_pointer_vectors{tt}(valid_reverse_vector_inds);
                % Find which filaments in the next frame should be forwarded to
                tracked_indices_array{tt+1}(valid_reverse_vector_inds) = ...
                    tracked_indices_array{tt}(valid_reverse_vector_pointers);
                
                % ----
                % Second, assign a new index to all filaments in the next frame
                % that are not connected to a filament in the last frame
                
                % Find invalid forwarding indices that need a new tracking
                % index
                newly_tracked_filaments = find(reverse_pointer_vectors{tt}==-1);
                %Number of filaments that need a new tracking index
                number_newly_tracked = numel(newly_tracked_filaments);
                if number_newly_tracked > 0
                    % Assign the filaments that need it new indices in the next
                    % frame
                    tracked_indices_array{tt+1}(newly_tracked_filaments) = ...
                        double(current_index_counter) + [1:number_newly_tracked];
                    % Increment the counter needed for increasing assignment of new
                    % indices
                    current_index_counter = ...
                        current_index_counter+number_newly_tracked;
                end
                
            end
            
            % Just a checksum to verify that the code isn't messed up
            if current_index_counter ~= ...
                    number_of_nonconnected_filaments + filaments_in_frame(1);
                sprintf(['Number of different tracking indices: %d\n,' ...
                    'Number of filament track counters set: %d'],...
                    current_index_counter, ...
                    number_of_nonconnected_filaments + filaments_in_frame(1))
                error('Something went wrong during index forwarding.')
            end
            number_of_tracked_filaments = current_index_counter;
            
            % ----
            % Create a control video of binary thresholding and filament
            % tracking
            if tracking_control_video
                
                % Make up a random color palette to use for video construction
                base_brightness = 0.4; % To avoid all black filaments
                below_top_brightness = 0.05; % To avoud all white filaments
                filament_colors = uint8( ...
                    255.*(...
                    (1-base_brightness).*ones(number_of_tracked_filaments,3) ...
                    + (base_brightness-below_top_brightness).*rand(number_of_tracked_filaments,3)));
                
                tracking_frames = zeros(height,width,1,number_of_merged_frames,'uint8');
                
                for ff = 1:number_of_merged_frames
                    this_frame_filaments = frame_conncomp_props{ff};
                    this_frame_image = 255.*ones(height,width,'uint8');
                    for nn = 1:filaments_in_frame(ff)
                        this_filament_pixel_list = ...
                            this_frame_filaments(nn).PixelIdxList;
                        color_from_tracking = filament_colors(...
                            tracked_indices_array{ff}(nn));
                        this_frame_image(this_filament_pixel_list) = ...
                            color_from_tracking;
                    end
                    tracking_frames(:,:,1,ff) = this_frame_image;
                end
                
                try %write the tracking control video
                    
                    % ----
                    % Open VideoWrite object
                    tracking_video = VideoWriter(...
                        [this_target_folder filesep 'Tracking_Control.avi'],...
                        'Uncompressed AVI');
                    tracking_video.FrameRate = frame_rate./frames_to_merge;
                    
                    % ----
                    %Open, write into, and close file associated with video object
                    open(tracking_video);
                    writeVideo(tracking_video,tracking_frames);
                    close(tracking_video)
                    % Remove tracking_video and tracking frames from memory
                    tracking_video = [];
                    tracking_frames = [];
                catch ME %make an entry into the error log
                    error_log{vv} = sprintf(...
                        'Error writing tracking video file :\n%s', ...
                        [this_target_folder filesep 'Tracking_Control.avi']);
                    no_error_flag(vv) = false;
                    fid = fopen([this_target_folder filesep 'error_report.txt'],'w');
                    fprintf(fid,'%s\n',error_log{vv});
                    fclose(fid);
                    continue %skip to next video in the batch
                end
                
            end
            
            % ----
            % Store filament properties in a structure array of size N-by-1,
            % where N is the number of tracked filaments
            
            % Stop analyzing this video and write an error report if no
            % filaments are left from this video
            if number_of_tracked_filaments == 0
                no_error_flag(vv) = false;
                error_log{vv} = sprintf(...
                    ['No traces could be detected, try relaxing analysis parameters, in\n',...
                    '\n%s'], ...,
                    [this_target_folder]);
                fid = fopen([this_target_folder filesep 'error_report.txt'],'w');
                fprintf(fid,'%s\n',error_log{vv});
                fclose(fid);
                continue
            end
            
            tracked_filaments = struct;
            tracked_filaments(number_of_tracked_filaments).appearance_time = [];
            tracked_filaments(end).appearance_frame = [];
            tracked_filaments(number_of_tracked_filaments).presence_time = [];
            tracked_filaments(end).present_frames = [];
            
            % Calculate for each tracked filament the time it is present for,
            % in units of seconds.
            
            indices_from_all_frames = [tracked_indices_array{:}];
            for nn = 1:number_of_tracked_filaments
                present_frames = sum(indices_from_all_frames==nn);
                tracked_filaments(nn).present_frames = present_frames;
                tracked_filaments(nn).presence_time = ...
                    present_frames.*frames_to_merge./frame_rate;
            end
            
            % ----
            % Use the tracked_indices_array to collect the properties of all
            % tracked filaments into the tracked_filaments structure array
            
            % To keep track how many frames have already been stored into a
            % specific filaments element of the structure array
            [tracked_filaments(:).frames_tracked] = ...
                deal(0);
            
            % Pre-allocate the fields that will contain per frame properties.
            % These fields will be numeric or cell arrays, and be of different
            % lengths, so they need to be allocated individually.
            
            tracked_filaments(end).filament_lengths = [];
            tracked_filaments(end).filament_widths = [];
            tracked_filaments(end).filament_centroid_xx = [];
            tracked_filaments(end).filament_centroid_yy = [];
            tracked_filaments(end).filament_pixel_indices_linear = [];
            tracked_filaments(end).filament_pixel_indices_xx = {};
            tracked_filaments(end).filament_pixel_indices_yy = {};
            tracked_filaments(end).filament_images = {};
            tracked_filaments(end).filament_solidities = [];
            tracked_filaments(end).filament_bounding_boxes = {};
            
            for nn = 1:number_of_tracked_filaments
                tracked_filaments(nn).filament_lengths = ...
                    zeros(1,tracked_filaments(nn).present_frames);
                tracked_filaments(nn).filament_widths = ...
                    zeros(1,tracked_filaments(nn).present_frames);
                tracked_filaments(nn).filament_centroid_xx = ...
                    zeros(1,tracked_filaments(nn).present_frames);
                tracked_filaments(nn).filament_centroid_yy = ...
                    zeros(1,tracked_filaments(nn).present_frames);
                tracked_filaments(nn).filament_pixel_indices_linear = ...
                    cell(1,tracked_filaments(nn).present_frames);
                tracked_filaments(nn).filament_pixel_indices_xx = ...
                    cell(1,tracked_filaments(nn).present_frames);
                tracked_filaments(nn).filament_pixel_indices_yy = ...
                    cell(1,tracked_filaments(nn).present_frames);
                tracked_filaments(nn).filament_images = ...
                    cell(1,tracked_filaments(nn).present_frames);
                tracked_filaments(nn).filament_solidities = ...
                    zeros(1,tracked_filaments(nn).present_frames);
                tracked_filaments(nn).filament_bounding_boxes = ...
                    cell(1,tracked_filaments(nn).present_frames);
            end
            
            % Assign the properties to the tracked filaments
            
            for tt = 1:number_of_merged_frames
                for ff = 1:filaments_in_frame(tt)
                    
                    % Find the index that the filament is tracked by
                    tracking_index = tracked_indices_array{tt}(ff);
                    % Increment the counter which keeps track of how many frames
                    % have already been stored for this specific tracked
                    % filament
                    
                    tracked_filaments(tracking_index).frames_tracked = ...
                        tracked_filaments(tracking_index).frames_tracked + 1;
                    % The "how-maniest" frame that the filament is tracked for
                    % this frame is
                    track_frame = ...
                        tracked_filaments(tracking_index).frames_tracked;
                    
                    % For the first frame that the filament is tracked for,
                    % store the frame number and the time of the frame
                    if tracked_filaments(tracking_index).frames_tracked == 1
                        tracked_filaments(tracking_index).appearance_frame = tt;
                        tracked_filaments(tracking_index).appearance_time = ...
                            tt.*frames_to_merge./frame_rate;
                    end
                    
                    % Assign properties
                    tracked_filaments(tracking_index).filament_lengths(track_frame) = ...
                        frame_conncomp_props{tt}(ff).Length;
                    tracked_filaments(tracking_index).filament_widths(track_frame) = ...
                        frame_conncomp_props{tt}(ff).Width;
                    tracked_filaments(tracking_index).filament_centroid_xx(track_frame) = ...
                        frame_conncomp_props{tt}(ff).WeightedCentroid(1);
                    tracked_filaments(tracking_index).filament_centroid_yy(track_frame) = ...
                        frame_conncomp_props{tt}(ff).WeightedCentroid(2);
                    tracked_filaments(tracking_index).filament_pixel_indices_linear{track_frame} = ...
                        frame_conncomp_props{tt}(ff).PixelIdxList;
                    tracked_filaments(tracking_index).filament_pixel_indices_xx{track_frame} = ...
                        frame_conncomp_props{tt}(ff).PixelList(:,1);
                    tracked_filaments(tracking_index).filament_pixel_indices_yy{track_frame} = ...
                        frame_conncomp_props{tt}(ff).PixelList(:,2);
                    tracked_filaments(tracking_index).filament_images{track_frame} = ...
                        frame_conncomp_props{tt}(ff).Image;
                    tracked_filaments(tracking_index).filament_solidities(track_frame) = ...
                        frame_conncomp_props{tt}(ff).Solidity;
                    tracked_filaments(tracking_index).filament_bounding_boxes{track_frame} = ...
                        frame_conncomp_props{tt}(ff).BoundingBox;
                end
            end
            
            % Check if all frames have been tracked properly
            tracked_filaments_check = [tracked_filaments(:).frames_tracked];
            tracked_indices_check = [tracked_filaments(:).present_frames];
            check_differences = tracked_filaments_check ...
                -  tracked_indices_check;
            if sum(check_differences)==0
                % Remove the field used to track how many frames have been added
                tracked_filaments = ...
                    rmfield(tracked_filaments,'frames_tracked');
            else
                error(...
                    'Something went wrong while assigning properties to tracked filaments.')
            end
            
            % Removal of filaments that have not been continuously followed for
            % a specified minimum presence time
            
            long_enough_presence_ids = ...
                find([tracked_filaments.presence_time]>=min_presence_time);
            tracked_filaments = tracked_filaments(long_enough_presence_ids);
            number_of_traces = numel(long_enough_presence_ids);
            
            % Stop analyzing this video and write an error report if no
            % filaments are left from this video
            if number_of_traces == 0
                no_error_flag(vv) = false;
                error_log{vv} = sprintf(...
                    ['No traces left after rejection, try relaxing parameters, in\n',...
                    '\n%s'], ...,
                    [this_target_folder]);
                fid = fopen([this_target_folder filesep 'error_report.txt'],'w');
                fprintf(fid,'%s\n',error_log{vv});
                fclose(fid);
                continue
            end
            
            %% Trace and frame-to-frame velocity, centroid trace construction
            
            for nn = 1:number_of_traces
                
                % Linear indices to the pixels that are part of this filament's
                % trace
                tracked_filaments(nn).trace_pixel_indices_linear = unique(...
                    cat(1,tracked_filaments(nn).filament_pixel_indices_linear{:}));
                
                % Binary image of filament trace. The image is of the size of a
                % whole frame, that is height-by-width
                trace_frame = false(height,width);
                true_inds = ind2sub([height,width], ...
                    tracked_filaments(nn).trace_pixel_indices_linear);
                trace_frame(true_inds) = true;
                tracked_filaments(nn).trace_frame = trace_frame;
                
                % Calculate region properties of the trace
                this_trace_properties = regionprops(trace_frame, ...
                    {'Area','Perimeter','Centroid', ...
                    'Image','Solidity','BoundingBox'});
                
                % Length and width calculation
                area = this_trace_properties.Area;
                perim = this_trace_properties.Perimeter;
                
                % Length and width of filament trace
                tracked_filaments(nn).trace_length = ...
                    rect_length(area,perim);
                tracked_filaments(nn).trace_width = ...
                    rect_width(area,perim);
                
                % Filament average length and width over all frames
                tracked_filaments(nn).average_filament_length = ...
                    mean([tracked_filaments(nn).filament_lengths(:)]);
                tracked_filaments(nn).average_filament_width = ...
                    mean([tracked_filaments(nn).filament_widths(:)]);
                
                % Trace velocity
                tracked_filaments(nn).trace_velocity = ...
                    ( tracked_filaments(nn).trace_length ...
                    - tracked_filaments(nn).average_filament_length) ...
                    ./tracked_filaments(nn).presence_time;
                %To avoid negative trace velocities due to subtraction of
                %average filament length
                tracked_filaments(nn).trace_velocity = ...
                    (tracked_filaments(nn).trace_velocity>=0) ...
                    .* tracked_filaments(nn).trace_velocity;
                
                
                % Centroid frame-to-frame velocities
                centroid_xx = tracked_filaments(nn).filament_centroid_xx;
                centroid_yy = tracked_filaments(nn).filament_centroid_yy;
                tracked_filaments(nn).frame_to_frame_velocities = ...
                    ( (centroid_xx(2:end)-centroid_xx(1:end-1)).^2 ...
                    + (centroid_yy(2:end)-centroid_yy(1:end-1)).^2).^0.5 ...
                    .*experimental_scaling_factor.*frame_rate./frames_to_merge;
                
                tracked_filaments(nn);
                
                % Store some of the region properties for this trace
                tracked_filaments(nn).trace_centroid_xx = ...
                    this_trace_properties.Centroid(1);
                tracked_filaments(nn).trace_centroid_yy = ...
                    this_trace_properties.Centroid(2);
                % Image of this trace inside of a tight fitting bounding box
                tracked_filaments(nn).trace_image = ...
                    this_trace_properties.Image;
                % Solidity of this trace
                tracked_filaments(nn).trace_solidity = ...
                    this_trace_properties.Solidity;
                % Bounding box of this trace
                tracked_filaments(nn).trace_bounding_box = ...
                    this_trace_properties.BoundingBox;
                
                
            end
            
            % Remove traces that are shorter than required minimum length
            long_enough_trace_ids = ...
                find(([tracked_filaments.trace_length] ...
                - [tracked_filaments.average_filament_length]) ...
                >= min_trace_length);
            tracked_filaments = tracked_filaments(long_enough_trace_ids);
            number_of_traces = numel(long_enough_trace_ids);
            
            % Stop analyzing this video and write an error report if no
            % filaments are left from this video
            if number_of_traces == 0
                no_error_flag(vv) = false;
                error_log{vv} = sprintf(...
                    ['No traces left after rejection, try relaxing parameters, in\n',...
                    '\n%s'], ...,
                    [this_target_folder]);
                fid = fopen([this_target_folder filesep 'error_report.txt'],'w');
                fprintf(fid,'%s\n',error_log{vv});
                fclose(fid);
                continue
            end
            
            
            
            %% Immotile frame removal (if switched on as an option)
                        
            if islogical(immotile_frame_removal) && immotile_frame_removal
                % If removal of frames below automatically determined
                % threshold velocity was selected
                
                try %Try immotile frame removal
                    
                    % Read in values that will be used to determine the threshold
                    % frame to frame velocity. This threshold will be used to
                    % detect immotile frames in filament motion
                    
                    imm_trace_velocities = ...
                        [tracked_filaments(:).trace_velocity];
                    imm_average_lengths = ...
                        [tracked_filaments(:).average_filament_length];
                    imm_frame_to_frame_velocities = ...
                        [tracked_filaments(:).frame_to_frame_velocities];
                    
                    % ----
                    % Determine a threshold that separates the motile from the
                    % immotile population
                    
                    % Histogram to determine region with sufficient filament counts
                    % to not confuse the curve fitting
                    velocities = real(imm_frame_to_frame_velocities);
                    hist_bin_number = 150;
                    minimum_count = 3;
                    hist_bin_edges = linspace(0,max(velocities),hist_bin_number+1);
                    hist_bin_centers = ...
                        (hist_bin_edges(1:end-1)+hist_bin_edges(2:end))./2;
                    [NN,bin_indices] = histc(velocities,hist_bin_edges);
                    min_count_inds = find(NN>minimum_count);
                    
                    % Collect velocities from the bins with sufficient count
                    velocities_trimmed = [];
                    for bb = min_count_inds
                        this_bin_inds = find(bin_indices==bb);
                        velocities_trimmed = ...
                            [velocities_trimmed,velocities(this_bin_inds)];
                    end
                    
                    % Fit two Gaussian mixture model
                    fit_options = statset('Display','off','MaxIter',2000, ...
                        'TolFun',1e-6);
                    mixture_fit = gmdistribution.fit( ...
                        velocities_trimmed.',2,'Options',fit_options);
                    mixture_model = @(vv) pdf(mixture_fit,vv);
                    vv_support = linspace(0,max(velocities_trimmed),300).';
                    frequencies = mixture_model(vv_support);
                    
                    % Two Gaussians mixture model parameters
                    component_proportions = mixture_fit.PComponents;
                    component_means = mixture_fit.mu;
                    component_sigma = squeeze(mixture_fit.Sigma).';
                    
                    % Sort parameters so that first component is that with lower mean
                    [component_means,sort_inds] = sort(component_means,'ascend');
                    component_proportions = component_proportions(sort_inds);
                    component_sigma = component_sigma(sort_inds);
                    motile_fraction = component_proportions(2);
                    
                    % Store mean of higher velocity Gaussian as v_max value
                    % from frame-to-frame velocities
                    
                    frame_to_frame_vmax = component_means(2);
                    
                    % Find minimum value between the two peaks
                    minimum_search_support = ...
                        linspace(component_means(1),component_means(2),2000).';
                    minimum_search_frequencies = ...
                        mixture_model(minimum_search_support);
                    [min_value,min_index] = min(minimum_search_frequencies);
                    %This is the actual threshold used for separation into motile
                    %and immotile
                    imm_threshold_velocity = minimum_search_support(min_index);
                    
                    
                    % Remove the immotile frames from the calculation of the
                    % trace velocity, remove traces that do not reach the
                    % minimum presence time after immotile frame removal
                    
                    for nn = 1:number_of_traces
                        % Number of all frames
                        all_present_frames = ...
                            tracked_filaments(nn).present_frames;
                        % Number of frames that are motile
                        motile_frames = sum( ...
                            tracked_filaments(nn).frame_to_frame_velocities ...
                            >imm_threshold_velocity);
                        % Reassign only motile frames as present frames
                        tracked_filaments(nn).present_frames = ...
                            motile_frames;
                        % Assign only motile time as presence time
                        tracked_filaments(nn).presence_time = ...
                            tracked_filaments(nn).present_frames ...
                            .*frames_to_merge./frame_rate;
                        % Adjust trace velocity
                        tracked_filaments(nn).trace_velocity = ...
                            tracked_filaments(nn).trace_velocity ...
                            .*all_present_frames./motile_frames;
                    end
                    
                    % Removal of filaments that have not been continuously followed for
                    % a specified minimum presence time, now also taking into
                    % account frames of immotility
                    
                    long_enough_presence_ids = ...
                        find([tracked_filaments.presence_time]>=min_presence_time);
                    tracked_filaments = tracked_filaments(long_enough_presence_ids);
                    number_of_traces = numel(long_enough_presence_ids);
                    
                    % Create a figure that will be saved to the target folder for
                    % control of the immotile frame removal procedure by the
                    % operator after analysis
                    
                    
                    % Open the figure
                    
                    current_figure = figure(vv);
                    
                    % Plot the trace velocities
                    
                    subplot(1,2,1)
                    % Trace velocities without removal
                    plot(real(imm_average_lengths),...
                        real(imm_trace_velocities),'ro')
                    hold on
                    % Line represting threshold velocity
                    plot([0 max(imm_average_lengths).*1.05], ...
                        imm_threshold_velocity.*ones(1,2),'k-')
                    %Trace velocities after removal of immotile frames
                    plot(real([tracked_filaments(:).average_filament_length]),...
                        real([tracked_filaments(:).trace_velocity]),'k+')
                    set(gca,'XLim',[0 max(imm_average_lengths).*1.05])
                    xlabel('Filament length[\mum]')
                    ylabel('Trace velocity[\mum/s]')
                    legend('Before frame removal',...
                        'Threshold velocity',...
                        'After frame removal')
                    hold off
                    
                    
                    subplot(1,2,2)
                    plot(hist_bin_centers,NN(1:end-1),'k-','LineWidth',3)
                    hold on
                    plot(hist_bin_centers(min_count_inds),NN(min_count_inds),...
                        'k--','LineWidth',1,'Color',[1 1 1])
                    plot(vv_support, ...
                        frequencies./max(frequencies).*max(NN),'r-')
                    plot(imm_threshold_velocity, ...
                        min_value./max(frequencies).*max(NN),'rs',...
                        'MarkerFaceColor',[1 0 0],'MarkerSize',8)
                    plot(imm_threshold_velocity.*ones(1,2),...
                        [0 max(NN).*1.05],'r-',...
                        'LineWidth',1.25)
                    set(gca,'XLim',[0 max(velocities_trimmed)])
                    set(gca,'YLim',[0 max(NN).*1.05])
                    ylabel('Frequency')
                    xlabel('Frame-to-frame velocity[\mum/s]')
                    legend('Frequency','Above min count','2 Gaussian model',...
                        'Threshold velocity','Threshold velocity')
                    hold off
                    
                    saveas(current_figure, [this_target_folder filesep 'Immotile_Frame_Removal_Check'], 'pdf')
                    
                    immotile_removal_success = true;
                    
                catch ME %Write an error to the error log batch
                    
                    disp(ME.getReport)
                    
                    immotile_removal_success = false;
                    error_log{vv} = sprintf(...
                        'Error removing immotile frames, affects target folder:\n%s\nResults are regular results without immotile frames removed.', ...
                        [this_target_folder]);
                    fid = fopen([this_target_folder filesep 'error_report.txt'],'w');
                    fprintf(fid,'%s\n',error_log{vv});
                    fclose(fid);
                    
                end
                
            elseif isnumeric(immotile_frame_removal) && ...
                    isfinite(immotile_frame_removal)
                
                % If removal of frames below manually set
                % threshold velocity was selected
                
                try %Try immotile frame removal
                    
                    % Read in values that will be used to determine the threshold
                    % frame to frame velocity. This threshold will be used to
                    % detect immotile frames in filament motion
                    
                    imm_trace_velocities = ...
                        [tracked_filaments(:).trace_velocity];
                    imm_average_lengths = ...
                        [tracked_filaments(:).average_filament_length];
                    
                    imm_threshold_velocity = immotile_frame_removal;
                    % Remove the immotile frames from the calculation of the
                    % trace velocity, remove traces that do not reach the
                    % minimum presence time after immotile frame removal
                    
                    for nn = 1:number_of_traces
                        % Number of all frames
                        all_present_frames = ...
                            tracked_filaments(nn).present_frames;
                        % Number of frames that are motile
                        motile_frames = sum( ...
                            tracked_filaments(nn).frame_to_frame_velocities ...
                            >imm_threshold_velocity);
                        % Reassign only motile frames as present frames
                        tracked_filaments(nn).present_frames = ...
                            motile_frames;
                        % Assign only motile time as presence time
                        tracked_filaments(nn).presence_time = ...
                            tracked_filaments(nn).present_frames ...
                            .*frames_to_merge./frame_rate;
                        % Adjust trace velocity
                        tracked_filaments(nn).trace_velocity = ...
                            tracked_filaments(nn).trace_velocity ...
                            .*all_present_frames./motile_frames;
                    end
                    
                    % Removal of filaments that have not been continuously followed for
                    % a specified minimum presence time, now also taking into
                    % account frames of immotility
                    
                    long_enough_presence_ids = ...
                        find([tracked_filaments.presence_time]>=min_presence_time);
                    tracked_filaments = tracked_filaments(long_enough_presence_ids);
                    number_of_traces = numel(long_enough_presence_ids);
                    
                    % Create a figure that will be saved to the target folder for
                    % control of the immotile frame removal procedure by the
                    % operator after analysis
                    
                    
                    % Open the figure
                    
                    current_figure = figure(vv);
                    
                    % Plot the trace velocities
                    
                    % Trace velocities without removal
                    plot(real(imm_average_lengths),real(imm_trace_velocities),'ro')
                    hold on
                    % Line represting threshold velocity
                    plot([0 max(imm_average_lengths).*1.05], ...
                        imm_threshold_velocity.*ones(1,2),'k-')
                    %Trace velocities after removal of immotile frames
                    plot(real([tracked_filaments(:).average_filament_length]),...
                        real([tracked_filaments(:).trace_velocity]),'k+')
                    set(gca,'XLim',[0 max(imm_average_lengths).*1.05])
                    xlabel('Filament length[\mum]')
                    ylabel('Trace velocity[\mum/s]')
                    legend('Before frame removal',...
                        'Threshold velocity',...
                        'After frame removal')
                    hold off
                    
                    saveas(current_figure, [this_target_folder filesep 'Immotile_Frame_Removal_Check'], 'pdf')
                    
                    immotile_removal_success = true;
                                        
                catch ME %Write an error to the error log batch
                    
                    immotile_removal_success = false;
                    error_log{vv} = sprintf(...
                        'Error removing immotile frames, affects target folder:\n%s\nResults are regular results without immotile frames removed.', ...
                        [this_target_folder]);
                    fid = fopen([this_target_folder filesep 'error_report.txt'],'w');
                    fprintf(fid,'%s\n',error_log{vv});
                    fclose(fid);
                    
                end
                
            end
            
            %% Save analysis results
            
            % ---
            
            %Store the breakage and trace results for this video into a cell
            %that contains data for the whole batch of videos that is analyzed
            %by this function run
            
            video_properties = struct;
            video_properties.width = width;
            video_properties.height = height;
            video_properties.duration = duration;
            video_properties.frame_rate = frame_rate;
            video_properties.merged_frame_rate = frame_rate./frames_to_merge;
            video_properties.traces = number_of_traces;
            video_properties.parameters = parameters(param_ind);
            video_properties.traces_checked = false;
            %Still add options here!!!
            
            if immotile_removal_success
                video_properties.frame_to_frame_vmax = frame_to_frame_vmax;
                video_properties.immotile_velocity_threshold = ...
                    imm_threshold_velocity;
                video_properties.motile_fraction = motile_fraction;
            end
            
            breakage_results = struct('video_properties',video_properties, ...
                'breakage_results',breakage_results);
            trace_results = struct('video_properties',video_properties, ...
                'trace_results',tracked_filaments);
            
            if trace_result_batch_requested
                trace_results_batch{vv} = trace_results;
            end
            if breakage_result_batch_requested
                breakage_results_batch{vv} = breakage_results;
            end
            
            % Send to save function, necessary for saving inside a
            % parallel for loop
            save_filepath = [this_target_folder filesep 'Analysis_Results.mat'];
            save_results(save_filepath,breakage_results,trace_results)
            
            % ----
            % Write to plain text file if requested by options
            
            if write_plain_text
                
                % Trace properties
                
                % Tab delimited file with all trace properties
                trace_text_file = ...
                    [this_target_folder filesep 'Trace_Results.dat'];
                names_of_fields = ...
                    {'trace_length',...
                    'trace_width',...
                    'trace_velocity',...
                    'trace_solidity'};
                field_headers = ...
                    {'Length[mum]',...
                    'Width[mum]',...
                    'Velocity[mum/s]',...
                    'Solidity'};
                
                number_of_columns = numel(names_of_fields);
                
                % Write the headers for data colums
                trace_result_file_writeobject = ...
                    fopen(trace_text_file,'w');
                for ff = 1:number_of_columns
                    fprintf(trace_result_file_writeobject,...
                        '%s\t',field_headers{ff});
                end
                fprintf(trace_result_file_writeobject,'\n');
                fclose(trace_result_file_writeobject);
                
                % Compile write array and write to file
                
                write_array = zeros(number_of_traces,number_of_columns);
                for ff = 1:number_of_columns
                    write_array(:,ff) = ...
                        real([tracked_filaments(:).(names_of_fields{ff})]);
                end
                dlmwrite(trace_text_file,write_array,'-append','delimiter','\t')
                
            end
            
        end

        fprintf('Video %d of %d analyzed\n',vv,videos_in_batch)
        
    catch ME
        
        %Executed if an error occurred during processing of the current
        %video. An error report is created based on the error message
        %issued due to the error that occurred during processing of this
        %video
        no_error_flag(vv) = false;
        error_log{vv} = sprintf(...
            'Error:%s, in\n%s\n', ...,
            ME.message,[this_target_folder]);
        fid = fopen([this_target_folder filesep 'error_report.txt'],'w');
        fprintf(fid,'%s\n',error_log{vv});
        fclose(fid);
        continue
        
    end
    
    
end

% Assign trace and breakage results to outputs
if trace_result_batch_requested
    varargout(1) = {trace_results_batch};
end
if breakage_result_batch_requested
    varargout(2) = {breakage_results_batch};
end

% List of folders saved to for filament dropper
if sum(no_error_flag)==0
    % In case all videos produced an error during analysis
    target_folder_list = {};
else
    % A list of paths to target folders in which the analysis ran fine
    % is produced as an output of the function
    target_folder_list = {videos{no_error_flag,2}};
end

% Error log preparation for output

%If no error occured return false as error log output
if sum(cellfun(@(xx) ~isempty(xx),error_log))==0
    error_log = false;
end

% Get elapsed time for this analysis procedure
elapsed_time = toc;
if elapsed_time_requested
    varargout{3} = elapsed_time;
end

function save_results(save_filepath,breakage_results,trace_results)

% MatLab save file with structure array containing all results
save(save_filepath,...
    'breakage_results','trace_results');