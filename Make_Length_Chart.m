% Pick a All_Results.mat file
[file,path] = uigetfile('Analysis_Results.mat','Pick an Analysis_Results.mat file');
full_path = [path file];

% Load the filament lengths, pixel indices, and video properties
trace_results = load(full_path,'trace_results');
trace_results = trace_results.trace_results;
video_properties = trace_results.video_properties;
trace_results = trace_results.trace_results;

filament_lengths = {trace_results.filament_lengths};
filament_centroids_xx = {trace_results.filament_centroid_xx};
filament_centroids_yy = {trace_results.filament_centroid_yy};
filament_image_pixels = {trace_results.filament_pixel_indices_linear};
number_of_filaments = numel(filament_lengths);


width = video_properties.width; height = video_properties.height;

filament_image = true(height,width);
for ff = 1:number_of_filaments
   filament_image(filament_image_pixels{ff}{1}) = false;
end
imagesc(filament_image);
colormap(gray)
for ff = 1:number_of_filaments
   text(filament_centroids_xx{ff}(1),filament_centroids_yy{ff}(1),...
       sprintf('%.2f',filament_lengths{ff}(1)))
end


title('Filament lengths in \mum')