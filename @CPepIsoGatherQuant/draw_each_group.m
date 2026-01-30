function draw_each_group(obj,raw_name,current_ratioMatrix,current_rts,...
    current_inten, low_mz_bound, high_mz_bound, selected_charge,...
    current_iso_rt_range, current_iso_name, dir_save, color_map, legend_map)
% Re-quantify each group
% input:
%   raw_name
%       the name of the raw (mgf) file
%   current_ratioMatrix
%       ratio matrix of quantification in current group
%   current_rts
%       retention time in current group
%   current_inten
%       intensity in current group
%   low_mz_bound
%       low precursor m/z bound
%   high_mz_bound
%       high precursor m/z bound
%   selected_charge
%       current precursor charge
%   current_iso_rt_range
%       retention times of current IMPs
%   current_iso_name
%       names of current IMPs
%   dir_save
%       the directory to save the plot
%   color_map
%       color map
%   legend_map
%       legend map


if nargin < 10
    legend_map = [];
end
if nargin < 9
    color_map = [];
end

rt_error_tol = 1; % RT match tolerance, choose 1 arbitrarily

% Sort MS1 signal (pair of retention time and intensity) by time
% Sort and denoise using a relative abundance threshold method
[sort_rts, sort_ratioMatrix, rt_grid, smoothed_intensity, ~, is_valid] = ...
    CQuantIMPGroupUtils.prepare_ms1_xic(...
        obj.m_cMs12DatasetIO, raw_name, current_rts, current_inten, current_ratioMatrix, ...
        obj.m_minMSMSnum, low_mz_bound, high_mz_bound, selected_charge);

if ~is_valid
    return;
end

% Extract the rt bound of XIC peak and convert to index bounds
[~, ~, is_skip_vec, peak_ranges] = ...
    CQuantIMPGroupUtils.prepare_peak_ranges_from_iso_rt_range(...
        rt_grid, current_iso_rt_range, rt_error_tol);

% Calculate the ratio on each XIC points using kernel method, and normalize
esti_ratio = CChromatogramUtils.calculate_kernel_ratio(rt_grid, sort_rts, sort_ratioMatrix, peak_ranges, false);

% Get deconvoluted XIC using revised RT
total_xic = {rt_grid, smoothed_intensity};
ric = CQuantIMPGroupUtils.build_ric_from_peaks(...
    rt_grid, smoothed_intensity, esti_ratio, peak_ranges, is_skip_vec);

% Check if the output directory exists
if ~exist(dir_save,'dir')
    mkdir(dir_save);
end

% Plot the total XIC and RIC
plot_xics(ric, current_iso_name, total_xic, dir_save, raw_name, ...
    low_mz_bound, high_mz_bound, selected_charge, color_map, legend_map);
end



function plot_xics(ric, current_iso_name, total_xic, dir_save, raw_name, ...
    low_mz_bound, high_mz_bound, selected_charge, color_map, legend_map)
% Plot the XIC of each IMP and the total XIC, all rt category
% input:
%   ric
%       retention times and intensities of each IMP
%   current_iso_name
%       names of current IMPs
%   total_xic
%       retention times, intensities of the total XIC
%   dir_save
%       the directory to save the plot
%   raw_name
%       the name of the raw (mgf) file
%   low_mz_bound
%       low precursor m/z bound
%   high_mz_bound
%       high precursor m/z bound
%   selected_charge
%       current precursor charge
%   color_map
%       color map
%   legend_map
%       legend map

% Sort the retention time intervals and categorize them
if any(cellfun(@(x) isempty(x), ric(:, 1)))
    del_rows = cellfun(@(x) isempty(x), ric(:, 1));
    ric(del_rows,:) = [];
    current_iso_name(del_rows,:) = [];
end
start_rt_array = cell2mat(cellfun(@(x) x(1), ric(:, 1), 'UniformOutput', false));
end_rt_array = cell2mat(cellfun(@(x) x(end), ric(:, 1), 'UniformOutput', false));
tolerance_rt = (end_rt_array-start_rt_array) ./ 5;
rt_intervals = [start_rt_array - tolerance_rt, end_rt_array + tolerance_rt];
[sorted_intervals, sort_idx] = sortrows(rt_intervals);
[categorized_intervals, categorized_indices] = categorize_intervals(sorted_intervals);
for idx_cat = 1:max(categorized_indices)
    % Extract the retention times and intensities of each IMP in the current category
    group_current_ric = ric(sort_idx(categorized_indices == idx_cat), :);
    group_current_iso_name = current_iso_name(sort_idx(categorized_indices == idx_cat));
    plot_each_xic_group(group_current_ric, total_xic, categorized_intervals{idx_cat}, group_current_iso_name, ...
        fullfile(dir_save,[raw_name, '_', num2str(low_mz_bound), '-', num2str(high_mz_bound), ...
        '_+', num2str(selected_charge), '_', num2str(idx_cat), '.svg']), color_map, legend_map);
end
end



function [categorized_intervals, categorized_indices] = categorize_intervals(intervals)
% Helper function to categorize retention time intervals
% input:
%   intervals
%       retention time intervals
% output:
%   categorized_intervals
%       categorized retention time intervals
%   categorized_indices
%       indices of the categories

categorized_intervals = {}; % Initialize an empty cell array to store categories
categorized_indices = zeros(length(intervals), 1); % Initialize an array to store the category index of each interval
for i = 1:size(intervals, 1)
    % Extract the current time interval
    currentInterval = intervals(i, :);
    
    % Initialize a flag to indicate if the current interval has been categorized
    categorized = false;
    
    % Iterate over the existing categories
    for j = 1:length(categorized_intervals)
        % Check if the current interval intersects with any interval in the category
        if is_intersecting(currentInterval, categorized_intervals{j})
            % If there is an intersection, add the current interval to this category
            categorized_intervals{j}(1) = min(currentInterval(1), categorized_intervals{j}(1)); %#ok<AGROW> 
            categorized_intervals{j}(2) = max(currentInterval(2), categorized_intervals{j}(2)); %#ok<AGROW> 
            categorized_indices(i) = j;
            categorized = true;
            break;
        end
    end
    
    % If the current interval does not intersect with any existing category, create a new category
    if ~categorized
        categorized_intervals{end+1} = currentInterval; %#ok<AGROW> 
        categorized_indices(i) = length(categorized_intervals);
    end
end
end



function intersection = is_intersecting(interval1, interval2)
% Helper function to check if two time intervals intersect
intersection = (interval1(1) <= interval2(2)) && (interval1(2) >= interval2(1));
end



function plot_each_xic_group(ric, total_xic, categorized_intervals, current_iso_name, file_name, color_map, legend_map)
% Plot the XIC of each IMP and the total XIC, only one rt category
% input:
%   ric
%       retention times and intensities of each IMP
%   total_xic
%       retention times, intensities of the total XIC
%   categorized_intervals
%       categorized retention time intervals [start_rt, end_rt]
%   current_iso_name
%       names of current IMPs
%   file_name
%       the file name to save the plot
%   color_map
%       color map
%   legend_map
%       legend map

% Init plot
% fig_w = 300;
fig_w = 2000;
fig_h = 800;
dpi = 300;
f = figure('Visible','off', 'Units','pixels', 'Position',[50, 50, fig_w, fig_h], 'Color','white');
set(f, 'PaperUnits','inches');
set(f, 'PaperPosition', [0, 0, fig_w/dpi, fig_h/dpi]);
set(f, 'PaperSize', [fig_w/dpi, fig_h/dpi]);
set(f, 'PaperPositionMode','manual', 'InvertHardcopy','off', 'Renderer','painters');
all_font_size = 7;
all_line_width = 1;
set(gca,'LooseInset',get(gca,'TightInset'), 'FontSize', all_font_size);
hold on

% Find the index of total boundary
total_idx_start = find(total_xic{1} <= categorized_intervals(1), 1, 'last');
total_idx_end = find(total_xic{1} >= categorized_intervals(2), 1);
if isempty(total_idx_start) || total_idx_start < 1
    total_idx_start = 1;
end
if isempty(total_idx_end) || total_idx_end > length(total_xic{1})
    total_idx_end = length(total_xic{1});
end
% Plot the total XIC
totalxic_plot = plot(total_xic{1}(total_idx_start:total_idx_end), total_xic{2}(total_idx_start:total_idx_end), ...
    'k', 'DisplayName','Total XIC');
set(totalxic_plot, 'LineWidth', all_line_width);

% Collect plot information for sorting
plot_info = struct('x_data', {}, 'y_data', {}, 'legend_string', {}, 'color', {});
plot_count = 0;

% Helper function to generate hash-based color from string to distinguish similar strings with different character orders
string_to_color = @(str) hsv2rgb([mod(hash_string_positional(str), 360)/360, 0.7, 0.9]);

% Collect information for each IMP plot
for idx_imp = 1:size(ric, 1)
    if trapz(ric{idx_imp, 1}, ric{idx_imp, 2}) < 1e-6
        continue;
    end
    
    plot_count = plot_count + 1;
    plot_info(plot_count).x_data = ric{idx_imp, 1};
    plot_info(plot_count).y_data = ric{idx_imp, 2};
    
    is_in_legend_map = ~isempty(legend_map) && legend_map.isKey(current_iso_name{idx_imp});
    is_in_color_map = ~isempty(color_map) && color_map.isKey(current_iso_name{idx_imp});

    if is_in_legend_map
        plot_info(plot_count).legend_string = ['XIC of ',legend_map(current_iso_name{idx_imp})];
    else
        legend_string = ['XIC of ',current_iso_name{idx_imp}];
        legend_string = strrep(legend_string, '_', '\_');
        legend_string = strrep(legend_string, '{', '\{');
        legend_string = strrep(legend_string, '}', '\}');
        plot_info(plot_count).legend_string = legend_string;
    end

    if is_in_color_map
        plot_info(plot_count).color = color_map(current_iso_name{idx_imp});
    else
        % Generate hash-based color from string
        plot_info(plot_count).color = string_to_color(current_iso_name{idx_imp});
    end
end

% Sort plot information by legend string
if plot_count > 0
    [~, sort_idx] = sort({plot_info.legend_string});
    plot_info = plot_info(sort_idx);
end

% Plot the XIC of each IMP in sorted order
for idx_plot = 1:plot_count
    isoxic_plot = plot(plot_info(idx_plot).x_data, plot_info(idx_plot).y_data, ...
        'DisplayName', plot_info(idx_plot).legend_string, 'Color', plot_info(idx_plot).color);
    set(isoxic_plot, 'LineWidth', all_line_width);
end

% Other setings of the plot
xlabel('Retention Time (min)', 'FontSize', all_font_size)
ylabel('Absolute intensity', 'FontSize', all_font_size)

ax = gca;
set(ax, 'Units','normalized');
set(ax, 'Position', [0.08 0.2 0.72 0.7]);
h_legend = legend('show', 'Location', 'northwest');
set(h_legend, 'FontSize', all_font_size, 'Box','off');
set(h_legend, 'Units','normalized');
set(h_legend, 'Position', [0.82 0.6 0.16 0.3]);
print(f, file_name, '-dsvg', ['-r', num2str(dpi)]);
end



% Position-sensitive string hash function
function hash_val = hash_string_positional(str)
    if isempty(str)
        hash_val = 0;
        return;
    end
    
    % Convert string to double array
    char_codes = double(str);
    
    % Method 1: Position-weighted polynomial hash (Horner's method)
    hash1 = 0;
    base = 31; % Prime number for better distribution
    for i = 1:length(char_codes)
        hash1 = mod(hash1 * base + char_codes(i), 2^32);
    end
    
    % Method 2: Position-indexed weighting
    positions = 1:length(char_codes);
    hash2 = sum(char_codes .* positions .* 137); % 137 is prime
    
    % Method 3: Alternating sign with position
    signs = (-1).^(positions - 1);
    hash3 = sum(char_codes .* signs .* positions);
    
    % Combine all three methods for maximum differentiation
    hash_val = mod(hash1 + hash2 * 17 + abs(hash3) * 23, 360);
end