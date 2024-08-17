% Shiying Li and Caroline Moosmueller, 2024. 
% plot SW relative error SW(sigma_k, mu)/SW(sigma_0, mu)to the target across different trials 
function plot_errors(sw_dist_all)

initialErr = sw_dist_all(1);
rel_sw_dist_all = sw_dist_all/initialErr; 

% Get the size of sw_dist_all
[num_trials, n1] = size(rel_sw_dist_all);

% Calculate average and standard deviation for full range
avg_error_full = mean(rel_sw_dist_all, 1); % Compute the average of each column (across trials)
std_error_full = std(rel_sw_dist_all, 0, 1); % Compute the standard deviation of each column (across trials)

% Plot individual curves (trials) - Full Range
figure; 

% Plot average curve with error bars for full range
xlabel('\boldmath$k$','Interpreter', 'latex', 'FontSize', 24, 'FontWeight', 'bold');

ylabel('\boldmath$\frac{SW_2(\sigma_k,\mu)}{SW_2(\sigma_0,\mu)}$', 'Interpreter', 'latex', 'FontSize', 25);


hold on; 
legend_entries = repmat({''}, 1, num_trials); % Preallocate cell array with empty strings

for i = 1:num_trials
    plot(0:n1-1, rel_sw_dist_all(i,:), '--','Color', [0.5, 0.5, 0.5], 'LineWidth', 1.5); 
end

errorbar(0:n1-1, avg_error_full, std_error_full, 'bdiamond-', 'LineWidth', 2.5, 'CapSize', 10, 'MarkerSize',6); % Plot average curve with error bars
legend_entries{num_trials} = 'Average \pm Std. Dev.'; % Add legend entry for average curve

h_legend = legend(legend_entries); % Create legend with dynamic entries
set(h_legend, 'FontSize', 16); % Set the font size of the legend
set(h_legend, 'Location', 'northwest'); % Set the location of the legend


set(gcf, 'Color', 'w')

hold off; % Disable hold to prevent further overlaying of plots


% Plot individual curves (trials) - Zoomed In
figure; 
hold on; 
initial_zoom_ind = 10; % for single is 10, for matrix is 2
% initial_zoom_ind = 2; % for single is 10, for matrix is 2
for i = 1:num_trials
    plot(initial_zoom_ind-1:n1-1, rel_sw_dist_all(i, initial_zoom_ind:end), '--','Color', [0.5, 0.5, 0.5], 'LineWidth', 1.5); % Plot each row from the second element to the end as a dotted line with stars at each data point
end

% Calculate average and standard deviation for zoomed-in plot
avg_error_zoomed = avg_error_full(initial_zoom_ind:end);
std_error_zoomed = std_error_full(initial_zoom_ind:end);


% Plot average curve with error bars for zoomed-in plot
errorbar(initial_zoom_ind-1:n1-1, avg_error_zoomed, std_error_zoomed, 'bdiamond-', 'LineWidth', 2.5, 'CapSize', 8);
hold off; 

set(gcf, 'Color', 'w')























