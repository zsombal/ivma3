function [L_crit_lo,L_crit_hi,slopes] = Length_Resolved_Pooled

clf

%[input_file,input_path] = uigetfile(pwd,'Select input bundle.','*.mat');

files = uipickfiles('FilterSpec',pwd);
number_of_bundles = numel(files);

L_min = 0;
L_max = 5;
windows = 30;
window_width = 0.115;
complex_in = true;
non_dropped = false;

% Analysis parameters
V_motile = 0.325;
frac_m_lower_threshold = 0.5;
frac_m_upper_threshold = 0.8;

% Containers to store length-resolved analysis results
cumulative_curves = cell(1,number_of_bundles);
mean_velocity = cell(1,number_of_bundles);
scaled_deviation = cell(1,number_of_bundles);
motile_fraction = cell(1,number_of_bundles);
L_crit_lo = zeros(1,number_of_bundles);
L_crit_hi = zeros(1,number_of_bundles);
slopes = zeros(number_of_bundles,4);
intercepts = zeros(number_of_bundles,4);

for bb = 1:number_of_bundles
    
    fprintf('Starting bundle %d of %d...',bb,number_of_bundles)
    
    
    % --- Make empirical cumulative probability density from frame-to-frame
    % velocities
    query_function = @(section) ...
        [section.trace_results.frame_to_frame_velocities];
    f2f_velocities = extract_by_keywords(files{bb}, ...
        [],[],[],query_function,non_dropped);
    [prob_density,support] = ecdf([f2f_velocities{:}]);
    cumulative_curves{bb} = cat(2,support,prob_density);
    
    
    % ---
    % Motile fraction
    
    filament_function = @(section) ...
        {section.trace_results.frame_to_frame_velocities};
    window_function = @(window_values) ...
        sum([window_values{:}]>=V_motile)./numel([window_values{:}]);
    
    [analysis_output,window_centers,window_counts] = apply_length_resolved( ...
        files{bb},[],L_min,L_max,windows,window_width,...
        complex_in,non_dropped, ...
        filament_function,window_function);
    analysis_output = [analysis_output{:}];
    motile_fraction{bb} = analysis_output;
    

        
    % ---
    % Detect where motile fraction crosses motile fraction lower and upper
    % threshold
    
    % Find indices of all windows with motile fraction above lower
    % threshold
    lower_thresh_ind = find(motile_fraction{bb}>=frac_m_lower_threshold);
    lower_thresh_ind = lower_thresh_ind(1);
    L_crit_lo(bb) = window_centers(lower_thresh_ind);
    % Find indices of all windows with motile fraction above lower
    % threshold
    upper_thresh_ind = find(motile_fraction{bb}>=frac_m_upper_threshold);
    upper_thresh_ind = upper_thresh_ind(1);
    L_crit_hi(bb) = window_centers(upper_thresh_ind);

    
        
    % ---
    % Mean trace velocity
    filament_function = @(section) ...
        [section.trace_results.trace_velocity];
    window_function = @(window_values) ...
        mean(window_values);
    
    [analysis_output,window_centers,window_counts] = apply_length_resolved( ...
        files{bb},[],L_min,L_max,windows,window_width,...
        complex_in,non_dropped, ...
        filament_function,window_function);
    mean_velocity{bb} = [analysis_output{:}];
    
    
    
    % ---
    % Scaled deviation
    
    filament_function = @scaled_deviation_filament_function;
    
    window_function = @(window) mean(cellfun( ...
        @(xx)(abs(xx(1)-mean(cellfun(@(xx)xx(1),window))).*sqrt(xx(2).*xx(3)))...
        ,window));
    
    [analysis_output,window_centers,window_counts] = apply_length_resolved( ...
        files{bb},[],L_min,L_max,windows,window_width,...
        complex_in,non_dropped, ...
        filament_function,window_function);
    analysis_output = [analysis_output{:}];
    scaled_deviation{bb} = analysis_output;
    
    fprintf('done.\n')
        
end


% ---
% Linear fit to determine the slope above L_crit_up for each bundle
for bb = 1:number_of_bundles
    
    linear_function = @(param,LL) param(2).*LL+param(1);
    initial_guess = [1 0];
    
    below_indices = find(window_centers<=L_crit_lo(bb));
    above_indices = find(window_centers>=L_crit_hi(bb));
    
    if numel(below_indices)>1
        % Mean velocity below critical length
        params = polyfit(window_centers(below_indices),...
            mean_velocity{bb}(below_indices),1);
        intercepts(bb,1) = params(2);
        slopes(bb,1) = params(1);
        % Scaled deviation below critical length
        params = ...
            polyfit(window_centers(below_indices),...
            scaled_deviation{bb}(below_indices),1);
        intercepts(bb,3) = params(2);
        slopes(bb,3) = params(1);
    end
    if numel(above_indices)>1
        % Mean velocity below critical length
        params = ...
            polyfit(window_centers(above_indices),...
            mean_velocity{bb}(above_indices),1);
        intercepts(bb,2) = params(2);
        slopes(bb,2) = params(1);
        % Scaled deviation below critical length
        params = ...
            polyfit(window_centers(above_indices),...
            scaled_deviation{bb}(above_indices),1);
        intercepts(bb,4) = params(2);
        slopes(bb,4) = params(1);
    end
    
end


% ---
% Plot results

% Frame-to-frame histograms
subplot(1,4,1)
for bb = 1:number_of_bundles
    stairs(cumulative_curves{bb}(:,1),cumulative_curves{bb}(:,2),...
        'k-','LineWidth',1.25)
    hold on
end
xlabel('V_{f2f}[\mum/s]')
ylabel('p(V_{f2f}<=V_{thr})')
title('Cumulative curve')
% set(gca,'YLim',[0 max([mean_velocity{:}]).*1.05])
% set(gca,'XLim',[window_centers(1) window_centers(end)])
hold off
    

% Mean trace velocity
subplot(1,4,2)
for bb = 1:number_of_bundles
    plot(window_centers,mean_velocity{bb},'k-')
    hold on
    plot(window_centers(window_centers<=L_crit_lo(bb)),...
        polyval([slopes(bb,1),intercepts(bb,1)],...
        window_centers(window_centers<=L_crit_lo(bb))),...
        'r--','LineWidth',1.25)
    plot(window_centers(window_centers>=L_crit_hi(bb)),...
        polyval([slopes(bb,2),intercepts(bb,2)],...
        window_centers(window_centers>=L_crit_hi(bb))),...
        'r--','LineWidth',1.25)
end
xlabel('L[\mum]')
ylabel('V[\mum/s]')
title('Mean trace velocity')
set(gca,'YLim',[0 max([mean_velocity{:}]).*1.05])
set(gca,'XLim',[window_centers(1) window_centers(end)])
hold off

% Scaled deviation
subplot(1,4,3)
for bb = 1:number_of_bundles
    plot(window_centers,scaled_deviation{bb},'k-')
    hold on
    plot(window_centers(window_centers<=L_crit_lo(bb)),...
        polyval([slopes(bb,3),intercepts(bb,4)],...
        window_centers(window_centers<=L_crit_lo(bb))),...
        'r--','LineWidth',1.25)
    plot(window_centers(window_centers>=L_crit_hi(bb)),...
        polyval([slopes(bb,4),intercepts(bb,4)],...
        window_centers(window_centers>=L_crit_hi(bb))),...
        'r--','LineWidth',1.25)
end
xlabel('L[\mum]')
ylabel('\deltav[\mum/s]')
title('Scaled deviation')
set(gca,'YLim',[0 max([scaled_deviation{:}]).*1.05])
set(gca,'XLim',[window_centers(1) window_centers(end)])
hold off

% Motile fraction
subplot(1,4,4)
for bb = 1:number_of_bundles
    plot(window_centers,motile_fraction{bb},'k-')
    hold on
end
xlabel('L[\mum]')
ylabel('f_m')
title('Motile fraction')
set(gca,'YLim',[0 max([motile_fraction{:}]).*1.05])
set(gca,'XLim',[window_centers(1) window_centers(end)])
hold off


function return_cell = ...
    scaled_deviation_filament_function(section)

trace_velocity = [section.trace_results.trace_velocity];
presence_time = [section.trace_results.presence_time];
average_filament_length = [section.trace_results.average_filament_length];

return_cell = ...
    arrayfun(@(velocity,time,length)[velocity,time,length],...
    trace_velocity,presence_time,average_filament_length,...
    'UniformOutput',false);