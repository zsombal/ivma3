function BlurVideo(videos,parameters)

%Check if the videos batch conforms with the requirements
video_batch_size = size(videos);
if numel(size(video_batch_size))~=2 || video_batch_size(2) ~= 2 || ~(iscell(videos))
    error('ivma3:EvaluateVideo:ImproperVideoBatch',...
        'The video batch must be 2-by-N cell array for N videos.')
end

%Number of videos in the batch
videos_in_batch = video_batch_size(1);

% Set up blurring filter
% Gaussian low pass filtering, hsize*hsize filter box used, sigma is
% standard deviation of Gaussian
blurring_filter = fspecial('gaussian', ...
    parameters.filtersize, parameters.filtersigma);

parfor vv = 1:videos_in_batch
    
    do_analysis = true;
        
    try %Try if video can be read
        current_video = VideoReader(videos{vv,1});
    catch ME %otherwise store error message and prevsnt further analysis
        do_analysis = false;
    end
    
    if do_analysis
        
        %% Creation of target directory to save to
        this_target_folder = videos{vv,2};
        [~,~,~] = mkdir(this_target_folder);
        
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
        blurred_frames = zeros(height,width,number_of_frames);
        
        % ---
        %Merge frames and get brightness medians
        for ff = 1:number_of_frames
            
            % -------------
            % Read video frames into matrix
            % Matrix containing all frame data for frames merged into this
            % merged frame, dimensions are:
            % image-height,image-width,RGB-channel,frame
            % Take mean of the three color channels to get gray scale frame
            % data, also scale from [0,255] to [0,1] range
            current_gray_frame = ...
                mean( ...
                read(current_video, ...
                [ff,ff]),3)./255
            
            blurred_frames(:,:,ff) = ...
                imfilter(current_gray_frame,blurring_filter);
            
        end
        
        % Remove current_video from memory
        current_video = [];
        
        % Control output of rescaled video)
        blurred_frames = reshape(blurred_frames,...
            height,width,1,number_of_frames);
        
        % ----
        % Open VideoWrite object
        blurred_video = VideoWriter(...
            [this_target_folder filesep 'BlurredVideo.avi'],...
            'Uncompressed AVI');
        blurred_video.FrameRate = frame_rate;
        
        % ----
        %Open, write into, and close file associated with video object
        open(blurred_video);
        %preprocesses_frames has to be scaled to [0,255] range and then
        %converted to unsigned integer to fit the uncompressed .avi
        %format for output
        writeVideo(blurred_video,uint8(blurred_frames.*255));
        close(blurred_video)
        % Remove preprocessing_video and preprocessed_frames
        % from memory
        blurred_video = [];
        blurred_frames = [];
                
    end
    
    fprintf('Video %d of %d blurred.\n',vv,videos_in_batch)
end