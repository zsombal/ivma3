% Algorithm parameters
false_pos_cost = 15; false_neg_cost = 1; % Costs of misclassification
ntrees = 150; % Number of decision trees in bag

% Open the bundles used for training and testing of the data classification
% system
[file_name, file_path] = uigetfile(pwd,...
    'Select bundle to train classification system.');
training_bundle = [file_path file_name];
[file_name, file_path] = uigetfile(file_path,...
    'Select bundle to test classification system.');
test_bundle = [file_path file_name];
[file_path file_name] = uiputfile(file_path,...
    'Select file to save trained model to.');
training_file = [file_name file_path];

non_dropped = true;
sufficient_keywords = {};
necessary_keywords = {};
not_keywords = {};

% ---
% Check if all videos have the same height and width, if not return error
% message

fprintf('Reading in video dimensions...')

% Get height and width of the different sections
get_height = @(section) ...
    section.video_properties.height;
get_width = @(section) ...
    section.video_properties.width;
training_heights = ...
    extract_by_keywords( training_bundle, ...
    sufficient_keywords,necessary_keywords,not_keywords, ...
    get_height,non_dropped);
training_widths = ...
    extract_by_keywords( training_bundle, ...
    sufficient_keywords,necessary_keywords,not_keywords, ...
    get_width,non_dropped);
test_heights =  ...
    extract_by_keywords( test_bundle, ...
    sufficient_keywords,necessary_keywords,not_keywords, ...
    get_height,non_dropped);
test_widths =  ...
    extract_by_keywords( test_bundle, ...
    sufficient_keywords,necessary_keywords,not_keywords, ...
    get_width,non_dropped);
heights = [training_heights{:} test_heights{:}];
widths = [training_widths{:} test_widths{:}];

fprintf('done.\n')

% Check if heights or widths differ between sections
if numel(unique(heights))>1 || numel(unique(widths))>1 %based on unique elements
    fprintf(['The videos have different dimensions.\n' ...
        'You can not train the system on videos of different dimensions.\n'])
    msgbox(['The videos have different dimensions.' ...
        'You can not train the system on videos of different dimensions.'],...
        'Video dimension mismatch','error')
else
    height = unique(heights);
    width = unique(widths);
end


% ---
% Extract tracedrop masks and trace image data and bounding boxes from
% results bundle

fprintf('Reading in training data...')

get_tracedrop_mask = @(section) ...
    section.tracedrop_mask;
get_trace_image = @(section) ...
    {section.trace_results.trace_image};
get_trace_bounding_box = @(section) ...
    {section.trace_results.trace_bounding_box};
get_trace_length = @(section) ...
    {section.trace_results.trace_length};
get_trace_width = @(section) ...
    {section.trace_results.trace_width};
get_trace_solidity = @(section) ...
    {section.trace_results.trace_solidity};

tracedrop_masks = ...
    extract_by_keywords( training_bundle, ...
    sufficient_keywords,necessary_keywords,not_keywords, ...
    get_tracedrop_mask,non_dropped);
trace_images = ...
    extract_by_keywords( training_bundle, ...
    sufficient_keywords,necessary_keywords,not_keywords, ...
    get_trace_image,non_dropped);
trace_bounding_boxes = ...
    extract_by_keywords( training_bundle, ...
    sufficient_keywords,necessary_keywords,not_keywords, ...
    get_trace_bounding_box,non_dropped);
trace_lengths = ...
    extract_by_keywords( training_bundle, ...
    sufficient_keywords,necessary_keywords,not_keywords, ...
    get_trace_length,non_dropped);
trace_widths = ...
    extract_by_keywords( training_bundle, ...
    sufficient_keywords,necessary_keywords,not_keywords, ...
    get_trace_width,non_dropped);
trace_solidities = ...
    extract_by_keywords( training_bundle, ...
    sufficient_keywords,necessary_keywords,not_keywords, ...
    get_trace_solidity,non_dropped);

tracedrop_masks = [tracedrop_masks{:}];
trace_images = [trace_images{:}];
trace_bounding_boxes = [trace_bounding_boxes{:}];
trace_lengths = [trace_lengths{:}];
trace_lengths = [trace_lengths{:}];
trace_widths = [trace_widths{:}];
trace_widths = [trace_widths{:}];
trace_solidities = [trace_solidities{:}];
trace_solidities = [trace_solidities{:}];

number_of_traces = numel(trace_images);

fprintf('done.\n%d traces in training data set.\n',number_of_traces)

%% Extract features for training data set

fprintf('Extracting trace features from training data set...')
colormap(gray)

% ---
%Specifications for corner detection
filter_coefficients = fspecial('gaussian',[21 1],2.5);
maximal_corners = 200;

% Determine number of featuresout-of-bundle data from
number_of_features = 2.*maximal_corners + 6;

% Extract and/or store trace features
feature_container = zeros(number_of_traces,number_of_features);
% Create trace images that are centered and have same 
centered_images = false(number_of_traces,height.*width);

parfor kk = 1:number_of_traces;
    
    
    % Initialize all 0 image
    this_image = false(height,width);
    
    % horizontal positioning
    xx_center = floor(trace_bounding_boxes{kk}(1) + ...
        trace_bounding_boxes{kk}(3)./2);
    xx_min = ceil(trace_bounding_boxes{kk}(1))-xx_center+round(width./2);
    xx_max = xx_min + trace_bounding_boxes{kk}(3)-1;
    % vertical positioning
    yy_center = floor(trace_bounding_boxes{kk}(2) + ...
        trace_bounding_boxes{kk}(4)./2);
    yy_min = ceil(trace_bounding_boxes{kk}(2))-yy_center+round(height./2);
    yy_max = yy_min + trace_bounding_boxes{kk}(4)-1;
    this_image(yy_min:yy_max,xx_min:xx_max) = ...
        trace_images{kk};
    
    % Get region properties of the trace in this image
    properties = regionprops(this_image,'Area','Orientation','Perimeter');
    
    % Rotate this image to have the trace main axis in horizontal direction
    this_image = imrotate(...
        this_image,-properties.Orientation,'bilinear','crop')
    
    %Store pixel data for this trace
    centered_images(kk,:) = this_image(:); %Store as row in centered_images
    
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
    
    feature_container(kk,:) = this_trace_features;
    
    % Plot the image of this trace with detected corners
    %     imagesc(this_image)
    %     hold on
    %     plot(corner_positions(:,1),corner_positions(:,2),'ro')
    %     hold off
    %     pause(1)
    
    
    
    fprintf('%d of %d traces done...\n',kk,number_of_traces)

end

feature_container(:,2.*maximal_corners+2) = trace_lengths;
feature_container(:,2.*maximal_corners+3) = trace_widths;
feature_container(:,2.*maximal_corners+4) = trace_solidities;

%Remove all imaginary parts from feature values
feature_container = real(feature_container);

% Clear data out of memory
clear trace_images trace_bounding_boxes heights widths

% ---
% Detect pixels for which any variation happens across the training data
% set
fprintf('Determining pixels containing any information...\n')

relevant_pixels = mean(centered_images,1);
relevant_pixels = relevant_pixels > 0 & relevant_pixels<1;

fprintf(...
    '%d of %d pixels contain any information useful for classification...\n', ...
    sum(relevant_pixels),numel(relevant_pixels))

cropped_pixel_data = double(centered_images(:,relevant_pixels));

% Clear data out of memory
clear centered_images

% Merge feature vector and relevant pixels vector
feature_container = cat(2,feature_container,cropped_pixel_data);
clear cropped_pixel_data

fprintf('done.\n')

% ---
% Train bagged decision tree method on training data set

% Create set of aggregated decision trees
fprintf('Creating bootstrap aggregation decision tree model...\n')
opt = statset('UseParallel','always');
decision_tree_model = ...
    TreeBagger(ntrees,feature_container,tracedrop_masks,...
    'oobpred','on','NPrint',1,'Options',opt);
fprintf('done.\n')

% ---
% Assess performance in training data set

fprintf('Assessing performance in training data set...\n')

% Display reduction of error in out-of-bag for increasing number of trees
figure(1)
subplot(1,3,1)
plot(oobError(decision_tree_model))
xlabel('Decision trees in ensemble')
ylabel('Prediction error rate')

% Make predictions inside the training data set
[predicted_class,prediction_score] = ...
    predict(decision_tree_model,feature_container);
predicted_class = logical(cellfun(@(xx) str2num(xx),predicted_class));

% True vs. false positive ROC curve
[false_pos,true_pos] = ...
    perfcurve(tracedrop_masks,prediction_score(:,2),true);
% subplot(1,3,2)
% plot(false_pos,true_pos,'LineWidth',1.25)
% title('ROC curve')out-of-bundle data from
% xlabel('False positive rate')
% ylabel('True positive rate')

% Cost optimization curve for false positives and false negatives
cost = [0 false_neg_cost;false_pos_cost 0];
[false_pos,e_cost,threshold] = ...
    perfcurve(tracedrop_masks,prediction_score(:,2),true,...
    'ycrit','ecost','cost',cost);
[min_cost,min_index] = min(e_cost);
subplot(1,3,2)
plot(threshold,e_cost./max(e_cost),...
    'k-','LineWidth',1.5)
hold on
plot(threshold(min_index),min_cost./max(e_cost),'ko')
plot([1 1].*threshold(min_index),[0 1],'k--')
plot(threshold,true_pos,'k-.','LineWidth',1.25)
plot(threshold,false_pos,'r--','LineWidth',1.25)
hold off
title('Cost optimization')
xlabel('Keep threshold')
ylabel('Estimated cost')

% Calculate and display confusion matrix
fprintf('Results in training data set:\n')
[training_confusion_matrix,class_order] = ...
    confusionmat(tracedrop_masks,predicted_class);
negative = training_confusion_matrix(1,1);
false_positive = training_confusion_matrix(1,2);
all_negative = negative+false_positive;
positive = training_confusion_matrix(2,2);
false_negative = training_confusion_matrix(2,1);
all_positive = positive + false_negative;
fprintf('%f3.3%% (%d of %d) dropped traces missed.\n',...
    100.*false_positive./all_negative,false_positive,all_negative)
fprintf('%f3.3%% (%d of %d) traces falsely dropped.\n',...
    100.*false_negative./all_positive,false_negative,all_positive)




% ---
% Test predictive power in second data set

% Extract tracedrop masks and trace image data and bounding boxes from
% results bundle

fprintf('Reading in test data...')

tracedrop_masks = ...
    extract_by_keywords( test_bundle, ...
    sufficient_keywords,necessary_keywords,not_keywords, ...
    get_tracedrop_mask,non_dropped);
trace_images = ...
    extract_by_keywords( test_bundle, ...
    sufficient_keywords,necessary_keywords,not_keywords, ...
    get_trace_image,non_dropped);
trace_bounding_boxes = ...
    extract_by_keywords( test_bundle, ...
    sufficient_keywords,necessary_keywords,not_keywords, ...
    get_trace_bounding_box,non_dropped);
trace_lengths = ...
    extract_by_keywords( test_bundle, ...
    sufficient_keywords,necessary_keywords,not_keywords, ...
    get_trace_length,non_dropped);
trace_widths = ...
    extract_by_keywords( test_bundle, ...
    sufficient_keywords,necessary_keywords,not_keywords, ...
    get_trace_width,non_dropped);
trace_solidities = ...
    extract_by_keywords( test_bundle, ...
    sufficient_keywords,necessary_keywords,not_keywords, ...
    get_trace_solidity,non_dropped);

tracedrop_masks = [tracedrop_masks{:}];
trace_images = [trace_images{:}];
trace_bounding_boxes = [trace_bounding_boxes{:}];
trace_lengths = [trace_lengths{:}];
trace_lengths = [trace_lengths{:}];
trace_widths = [trace_widths{:}];
trace_widths = [trace_widths{:}];
trace_solidities = [trace_solidities{:}];
trace_solidities = [trace_solidities{:}];

number_of_traces = numel(trace_images);

%% Extract features for test data set

fprintf('Extracting trace features from test data set...')
colormap(gray)

% Extract and/or store trace features
feature_container = zeros(number_of_traces,number_of_features);
% Create trace images that are centered and have same 
centered_images = false(number_of_traces,height.*width);

parfor kk = 1:number_of_traces;
    
    
    % Initialize all 0 image
    this_image = false(height,width);
    
    % horizontal positioning
    xx_center = floor(trace_bounding_boxes{kk}(1) + ...
        trace_bounding_boxes{kk}(3)./2);
    xx_min = ceil(trace_bounding_boxes{kk}(1))-xx_center+round(width./2);
    xx_max = xx_min + trace_bounding_boxes{kk}(3)-1;
    % vertical positioning
    yy_center = floor(trace_bounding_boxes{kk}(2) + ...
        trace_bounding_boxes{kk}(4)./2);
    yy_min = ceil(trace_bounding_boxes{kk}(2))-yy_center+round(height./2);
    yy_max = yy_min + trace_bounding_boxes{kk}(4)-1;
    this_image(yy_min:yy_max,xx_min:xx_max) = ...
        trace_images{kk};
    
    % Get region properties of the trace in this image
    properties = regionprops(this_image,'Area','Orientation','Perimeter');
    
    % Rotate this image to have the trace main axis in horizontal direction
    this_image = imrotate(...
        this_image,-properties.Orientation,'bilinear','crop')
    
    %Store pixel data for this trace
    centered_images(kk,:) = this_image(:); %Store as row in centered_images
    
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
    
    feature_container(kk,:) = this_trace_features;
    
    % Plot the image of this trace with detected corners
    %     imagesc(this_image)
    %     hold on
    %     plot(corner_positions(:,1),corner_positions(:,2),'ro')
    %     hold off
    %     pause(1)
    
    
    
    fprintf('%d of %d traces done...\n',kk,number_of_traces)

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
[predicted_class,prediction_score] = predict(decision_tree_model,feature_container);
predicted_class = logical(cellfun(@(xx) str2num(xx),predicted_class));
fprintf('done.\n')

%Calculate and display confusion matrix
fprintf('Results in test data set:\n')
test_confusion_matrix = confusionmat(tracedrop_masks,predicted_class);
negative = test_confusion_matrix(1,1);
false_positive = test_confusion_matrix(1,2);
all_negative = negative+false_positive;
positive = test_confusion_matrix(2,2);
false_negative = test_confusion_matrix(2,1);
all_positive = positive + false_negative;
fprintf('%f3.3%% (%d of %d) dropped traces missed.\n',...
    100.*false_positive./all_negative,false_positive,all_negative)
fprintf('%f3.3%% (%d of %d) traces falsely dropped.\n',...
    100.*false_negative./all_positive,false_negative,all_positive)


% True vs. false positive ROC curve
[false_pos,true_pos] = ...
    perfcurve(tracedrop_masks,prediction_score(:,2),true);
% subplot(1,3,3)
% plot(false_pos,true_pos,'LineWidth',1.25)
% title('ROC curve')
% xlabel('False positive rate')
% ylabel('True positive rate')

% Cost optimization curve for false positives and false negatives
[false_pos,e_cost,threshold] = ...
    perfcurve(tracedrop_masks,prediction_score(:,2),true,...
    'ycrit','ecost','cost',cost);
[min_cost,min_index] = min(e_cost);
opt_thresh = threshold(min_index);
subplot(1,3,3)
plot(threshold,e_cost./max(e_cost),...
    'k-','LineWidth',1.5)
hold on
plot(threshold(min_index),min_cost./max(e_cost),'ko')
plot([1 1].*opt_thresh,[0 1],'k--')
plot(threshold,true_pos,'k-.','LineWidth',1.25)
plot(threshold,false_pos,'r--','LineWidth',1.25)
hold off
title('Cost optimization')
xlabel('Keep threshold')
ylabel('Estimated cost')

%% Save model training file

save(training_file,'decision_tree_model','relevant_pixels','opt_thresh',...
    'height','width','-v7.3')

% ---
% Make a (V,L) plot displaying the manual and predicted data points
figure(2)

% Predicted mask
thresh_class_predictors = prediction_score(:,2)>=opt_thresh;
fprintf('Reading in (V,L) data for training data set...')

get_average_length = @(section) ...
    [section.trace_results.average_filament_length];
get_trace_velocity = @(section) ...
    [section.trace_results.trace_velocity];

average_lengths = ...
    extract_by_keywords( test_bundle, ...
    sufficient_keywords,necessary_keywords,not_keywords, ...
    get_average_length,non_dropped);
trace_velocities = ...
    extract_by_keywords( test_bundle, ...
    sufficient_keywords,necessary_keywords,not_keywords, ...
    get_trace_velocity,non_dropped);

average_lengths = [average_lengths{:}];
trace_velocities = [trace_velocities{:}];

fprintf('done.\n')

% (V,L diagram)
subplot(2,2,1)
plot(average_lengths(tracedrop_masks),...
    trace_velocities(tracedrop_masks),'ks')
hold on
plot(average_lengths(thresh_class_predictors),...
    trace_velocities(thresh_class_predictors),'r+')
hold off
xlabel('L[\mum]'),ylabel('V[\mum/s]')
title('All')

% 2D histogram difference plot
subplot(2,2,2)
v_max = max(trace_velocities);
L_max = max(average_lengths);
Edges = {linspace(0,L_max,25),linspace(0,v_max,10)};
[manual_counts,positions] = ...
    hist3(real([average_lengths(tracedrop_masks).',...
    trace_velocities(tracedrop_masks).']),...
    'Edges',Edges);
[predicted_counts,positions] = ...
    hist3(real([average_lengths(thresh_class_predictors).',...
    trace_velocities(thresh_class_predictors).']),...
    'Edges',Edges);
imagesc(positions{1},positions{2},...
    (manual_counts-predicted_counts)./manual_counts)
set(gca,'YDir','normal')
colorbar
xlabel('L[\mum]')
ylabel('V[\mum/s]')
title('(Pred-Man)/Man, all')


% (V,L) diagram and 2D difference histogram for only strictly real results
manual_average_lengths = average_lengths(tracedrop_masks);
manual_trace_velocities = trace_velocities(tracedrop_masks);
manual_real_indices = find(...
    ~imag(manual_average_lengths) & ...
    ~imag(manual_trace_velocities));
predicted_average_lengths = average_lengths(thresh_class_predictors);
predicted_trace_velocities = trace_velocities(thresh_class_predictors);
predicted_real_indices = find(...
    ~imag(predicted_average_lengths) & ...
    ~imag(predicted_trace_velocities));

% (V,L diagram)
subplot(2,2,3)
plot(manual_average_lengths(manual_real_indices),...
    manual_trace_velocities(manual_real_indices),'ks')
hold on
plot(predicted_average_lengths(predicted_real_indices),...
    predicted_trace_velocities(predicted_real_indices),'r+')
hold off
xlabel('L[\mum]'),ylabel('V[\mum/s]')
title('Only real')

% 2D histogram difference plot
subplot(2,2,4)
v_max = max(trace_velocities);
L_max = max(average_lengths);
Edges = {linspace(0,L_max,25),linspace(0,v_max,10)};
[real_manual_counts,real_positions] = ...
    hist3([manual_average_lengths(manual_real_indices).',...
    manual_trace_velocities(manual_real_indices).'],...
    'Edges',Edges);
[real_predicted_counts,real_positions] = ...
    hist3([predicted_average_lengths(predicted_real_indices).',...
    predicted_trace_velocities(predicted_real_indices).'],...
    'Edges',Edges);
imagesc(real_positions{1},real_positions{2},...
    (real_manual_counts-real_predicted_counts)./real_manual_counts)
set(gca,'YDir','normal')
colorbar
xlabel('L[\mum]')
ylabel('V[\mum/s]')
title('(Pred-Man)/Man, only real')


%% Display the false positive traces for user if desired
if strcmp(questdlg('Do you want to see the bad tracs that were missed?'),...
        'Yes')
    fprintf('Loading trace images for display...')
    trace_images = ...
        extract_by_keywords( test_bundle, ...
        sufficient_keywords,necessary_keywords,not_keywords, ...
        get_trace_image,non_dropped);
    trace_images = [trace_images{:}];
    trace_bounding_boxes = ...
        extract_by_keywords( test_bundle, ...
        sufficient_keywords,necessary_keywords,not_keywords, ...
        get_trace_bounding_box,non_dropped);
    trace_bounding_boxes = [trace_bounding_boxes{:}];
    fprintf('done.\n')
    false_positive_indices = find(thresh_class_predictors==true&...
        tracedrop_masks.'==false);
    figure(3)
    clf
    for kk = false_positive_indices.';
        % Initialize all 0 image
        this_image = false(height,width);
        
        % horizontal positioning
        xx_center = floor(trace_bounding_boxes{kk}(1) + ...
            trace_bounding_boxes{kk}(3)./2);
        xx_min = ceil(trace_bounding_boxes{kk}(1))-xx_center+round(width./2);
        xx_max = xx_min + trace_bounding_boxes{kk}(3)-1;
        % vertical positioning
        yy_center = floor(trace_bounding_boxes{kk}(2) + ...
            trace_bounding_boxes{kk}(4)./2);
        yy_min = ceil(trace_bounding_boxes{kk}(2))-yy_center+round(height./2);
        yy_max = yy_min + trace_bounding_boxes{kk}(4)-1;
        this_image(yy_min:yy_max,xx_min:xx_max) = ...
            trace_images{kk};
        colormap(gray)
        imagesc(this_image)
        pause(1)
    end
end
