function Plot_Length_Resolved

clf

keywords = [];

[input_file,input_path] = uigetfile(pwd,'Select input bundle.','*.mat');
input_bundle = [input_path input_file];
L_min = 0;
L_max = 10;
windows = 20;
window_width = 0.125;
complex_in = true;
non_dropped = false;

% ---
% Regular (V,L) diagram

%% Get the average filament lengths and the trace velocities

% Pull trace velocities out of the queried result sections
pull_velocity = ...
    @(result_section) [result_section.trace_results.trace_velocity];
trace_velocity = extract_by_keywords(input_bundle, ...
    [],[],[], ...
    pull_velocity,non_dropped);
trace_velocity = [trace_velocity{:}];

% Pull filament average lenghts out of the queried result sections
pull_average_length = ...
    @(result_section) ...
    [result_section.trace_results.average_filament_length];
average_filament_length = extract_by_keywords(input_bundle, ...
    [],[],[], ...
    pull_average_length,non_dropped);
average_filament_length = [average_filament_length{:}];

inds = 1:numel(average_filament_length);

if ~exist('non_dropped','var') || ~isempty(non_dropped) || ...
        ~non_dropeed
    non_dropped = false;
end

if ~complex_in
    % Keep only indices with real valued results in filament length and
    % filament trace velocity
   inds = find(imag(average_filament_length)==0 ...
    & imag(trace_velocity)==0);
end

subplot(2,3,1)
plot(average_filament_length(inds),trace_velocity(inds),'ko',...
    'MarkerSize',4)
xlabel('L[\mum]')
ylabel('V_{f2f}[\mum/s]')
title('Trace velocity')

% ---
% 2D histogram of f2f velocities and number of filaments

V_high = 1.5;
bins = 25;
binning_edges = linspace(0,V_high,bins+1);
bin_centers = (binning_edges(1:end-1)+binning_edges(2:end))./2;

% Filament and box car window level extraction functions
filament_function = @(section) ...
    {section.trace_results.frame_to_frame_velocities};
window_function = @(window_values) ...
    f2f_window_histogram(window_values,V_high,bins);

[analysis_output,window_centers,window_counts] = apply_length_resolved( ...
    input_bundle,[],L_min,L_max,windows,window_width, ...
    complex_in,non_dropped, ...
    filament_function,window_function);
normalized_frequencies = [analysis_output{:}];

subplot(2,3,2)
imagesc(window_centers,bin_centers,normalized_frequencies)
colormap(gray)
set(gca,'YDir','normal')
xlabel('L[\mum]')
title('Frame-to-frame velocity')
ylabel('V[\mum/s]')

subplot(2,3,3)
plot(window_centers,window_counts);
xlabel('L[\mum]')
ylabel('N filaments')
title('Filament count')
set(gca,'YLim',[0 max(window_counts).*1.05])
set(gca,'XLim',[0 L_max.*1.05])


% ---
% Mean trace velocity
filament_function = @(section) ...
    [section.trace_results.trace_velocity];
window_function = @(window_values) ...
    mean(window_values);

[analysis_output,window_centers,window_counts] = apply_length_resolved( ...
    input_bundle,[],L_min,L_max,windows,window_width,...
    complex_in,non_dropped, ...
    filament_function,window_function);
mean_trace_velocity = [analysis_output{:}];

subplot(2,3,4)
plot(window_centers,mean_trace_velocity);
xlabel('L[\mum]')
ylabel('V[\mum/s]')
title('Mean trace velocity')
set(gca,'YLim',[0 max(mean_trace_velocity).*1.05])
set(gca,'XLim',[0 L_max.*1.05])



% ---
% Scaled deviation

filament_function = @scaled_deviation_filament_function;
    
window_function = @(window) mean(cellfun( ...
    @(xx)(abs(xx(1)-mean(cellfun(@(xx)xx(1),window))).*sqrt(xx(2).*xx(3)))...
    ./(mean(cellfun(@(xx)xx(1),window))),window));
    
[analysis_output,window_centers,window_counts] = apply_length_resolved( ...
    input_bundle,[],L_min,L_max,windows,window_width,...
    complex_in,non_dropped, ...
    filament_function,window_function);
analysis_output = analysis_output;
scaled_deviation = [analysis_output{:}];

subplot(2,3,5)
plot(window_centers,scaled_deviation);
xlabel('L[\mum]')
ylabel('\deltav[\mum/s]')
title('Scaled deviation')
set(gca,'YLim',[0 max(scaled_deviation).*1.05])
set(gca,'XLim',[0 L_max.*1.05])


% ---
% Frequency power spectrum using Welch's method

%Get mean frame rate from results bundle
% Draw figure with spectra of all traces
pull_frame_rate = @(section) section.video_properties.merged_frame_rate;
fprintf(['Frame rates of different sections,\n' ...
    'Bad if they are different!\n'])
frame_rate = extract_by_keywords(input_bundle, ...
    [],[],[], ...
    pull_frame_rate,non_dropped);
frame_rate = [frame_rate{:}];
disp(frame_rate)
frame_rate = mean(frame_rate(:));

% Parameters for determing the power spectrum
frequency_resolution = 10;
welch_window_width = 2.*(frequency_resolution-1);
welch_window_overlap = []; fft_points = welch_window_width;
sampling_frequency = frame_rate;

% Filament and box car window level extraction functions
filament_function = @(section) ...
    {section.trace_results.frame_to_frame_velocities};
window_function = @(window_values) power_spectrum(window_values,...
    welch_window_width,welch_window_overlap,...
    fft_points,sampling_frequency);

% Call to length-resolved extraction function
[analysis_output,window_centers,window_counts] = apply_length_resolved( ...
    input_bundle,[],L_min,L_max,windows,window_width, ...
    complex_in,non_dropped, ...
    filament_function,window_function);

% Find box car windows for which power spectrum could be determined
successful_windows = find(cellfun(@(window) ...
    window.success==true,analysis_output));

%Plot the successfull spectra
subplot(2,3,6)

spectrum_grid = zeros(numel(successful_windows),frequency_resolution);
successful_centers = window_centers(successful_windows);
frequencies = analysis_output{successful_windows(1)}.frequencies;
[f_grid,L_grid] = meshgrid(frequencies,successful_centers);

for kk = 1:numel(successful_windows)
    spectrum_grid(kk,:) = ...
        analysis_output{successful_windows(kk)}.power_spectrum;
end

surfc(L_grid,f_grid,spectrum_grid)
xlabel('L[\mum]')
ylabel('f[1/s]')
zlabel('I[dB]')
axis tight
title('Frequency power spectrum')
hold off



function normalized_counts = ...
    f2f_window_histogram(window_values,V_high,bins)
% Gets the histogram of frame-to-frame velocities inside the window

binning_edges = linspace(0,V_high,bins+1);

window_values = [window_values{:}];

normalized_counts = histc(window_values,binning_edges);
normalized_counts = normalized_counts(1:end-1);
normalized_counts = normalized_counts./sum(normalized_counts);
normalized_counts = normalized_counts.';

function return_struct = ...
    power_spectrum(window_values,window_width,window_overlap,...
    fft_points,sampling_frequency)


power_spectrum_array = [];

per_frame_velocities = window_values;

long_enough_velocity_traces = 0;
for kk = 1:numel(per_frame_velocities)
    this_velocities = ...
        real(per_frame_velocities{kk});
    if numel(this_velocities)>window_width
        long_enough_velocity_traces = long_enough_velocity_traces+1;
        [power_spectrum,frequencies] = ...
            pwelch(this_velocities,window_width, ...
            window_overlap,fft_points,sampling_frequency);
        power_spectrum_array = [power_spectrum_array;power_spectrum.'];
    end
end

if long_enough_velocity_traces > 0
    % If at least one filament could be analyzed for its frequency spectrum
    
    power_spectrum_mean = mean(power_spectrum_array,1);
    power_to_dezibel = @(power_value) 10.*log10(power_value);
    power_spectrum_mean = power_to_dezibel(power_spectrum_mean);
    return_struct = struct;
    return_struct.power_spectrum = power_spectrum_mean;
    return_struct.frequencies = frequencies;
    return_struct.success = true;
else
    return_struct.success = false;
end


function return_cell = ...
    scaled_deviation_filament_function(section)

trace_velocity = [section.trace_results.trace_velocity];
presence_time = [section.trace_results.presence_time];
average_filament_length = [section.trace_results.average_filament_length];

return_cell = ...
    arrayfun(@(velocity,time,length)[velocity,time,length],...
    trace_velocity,presence_time,average_filament_length,...
    'UniformOutput',false);