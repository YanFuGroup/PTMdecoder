function drawGroup(ms12DatasetIO, minMSMSnum, raw_name, ratio_raw, rt_raw, ...
    intensity_raw, low_mz_bound, high_mz_bound, selected_charge, ...
    current_imp_rt_range, current_imp_name, dir_save, color_map, legend_map)
% Re-quantify each group
% input:
%   ms12DatasetIO (CMs12DatasetIO)
%       the dataset IO for accessing MS1 and MS2 spectra
%   minMSMSnum (1 x 1 double/int)
%       minimum number of MSMS spectra for a peptide to be considered
%   raw_name (1 x 1 char/string)
%       the name of the raw (mgf) file
%   ratio_raw (N x K double)
%       ratio matrix of quantification in current group
%   rt_raw (N x 1 double) minutes
%       retention time in current group
%   intensity_raw (N x 1 double) intensity
%       intensity in current group
%   low_mz_bound (1 x 1 double) m/z
%       low precursor m/z bound
%   high_mz_bound (1 x 1 double) m/z
%       high precursor m/z bound
%   selected_charge (1 x 1 double/int)
%       current precursor charge
%   current_imp_rt_range (K x 1 cell)
%       RT ranges for each IMP (each cell: [] or [rt_start, rt_end] in minutes)
%   current_imp_name (K x 1 cellstr/string)
%       names of current IMPs
%   dir_save (1 x 1 char/string)
%       the directory to save the plot
%   color_map (containers.Map or [])
%       color map (key: imp name, value: RGB 1x3)
%   legend_map (containers.Map or [])
%       legend map (key: imp name, value: display string)
    if nargin < 14
        legend_map = [];
    end
    if nargin < 13
        color_map = [];
    end

    rt_error_tol = 1; % RT match tolerance

    % Sort and denoise using a relative abundance threshold method
    [rt_sorted, ratio_sorted, xic_rt, xic_intensity_smoothed, ~, is_valid] = ...
        CXICPreprocessUtils.prepare_ms1_xic(...
            ms12DatasetIO, raw_name, rt_raw, intensity_raw, ratio_raw, ...
            minMSMSnum, low_mz_bound, high_mz_bound, selected_charge);

    if ~is_valid
        return;
    end

    % Extract the rt bound of XIC peak and convert to index bounds
    [~, is_skip_vec, xic_peak_idx_bounds] = ...
        CXICPreprocessUtils.prepare_peak_ranges_from_imp_rt_range(...
            xic_rt, current_imp_rt_range, rt_error_tol);

    % Calculate the ratio on each XIC point using kernel method, and normalize
    xic_ratio_estimated = CXICPeakUtils.calculate_kernel_ratio(...
        xic_rt, rt_sorted, ratio_sorted, xic_peak_idx_bounds, false);

    % Get deconvoluted XIC using revised RT
    total_xic = {xic_rt, xic_intensity_smoothed};
    ric = CXICAreaUtils.build_ric_from_peaks(...
        xic_rt, xic_intensity_smoothed, xic_ratio_estimated, ...
        xic_peak_idx_bounds, is_skip_vec);

    if ~exist(dir_save, 'dir')
        mkdir(dir_save);
    end

    plot_xics(ric, current_imp_name, total_xic, dir_save, raw_name, ...
        low_mz_bound, high_mz_bound, selected_charge, color_map, legend_map);
end



function plot_xics(ric, current_imp_name, total_xic, dir_save, raw_name, ...
    low_mz_bound, high_mz_bound, selected_charge, color_map, legend_map)
% Plot the XIC of each IMP and the total XIC, all rt category
% input:
%   ric (K x 2 cell)
%       retention times and intensities of each IMP; ric{i,1}=rt (minutes), ric{i,2}=intensity
%   current_imp_name (K x 1 cellstr/string)
%       names of current IMPs
%   total_xic (1 x 2 cell)
%       total XIC; total_xic{1}=rt (minutes), total_xic{2}=intensity
%   dir_save (1 x 1 char/string)
%       the directory to save the plot
%   raw_name (1 x 1 char/string)
%       the name of the raw (mgf) file
%   low_mz_bound (1 x 1 double) m/z
%       low precursor m/z bound
%   high_mz_bound (1 x 1 double) m/z
%       high precursor m/z bound
%   selected_charge (1 x 1 double/int)
%       current precursor charge
%   color_map (containers.Map or [])
%       color map
%   legend_map (containers.Map or [])
%       legend map

if any(cellfun(@(x) isempty(x), ric(:, 1)))
    del_rows = cellfun(@(x) isempty(x), ric(:, 1));
    ric(del_rows, :) = [];
    current_imp_name(del_rows, :) = [];
end
rt_start_array = cell2mat(cellfun(@(x) x(1), ric(:, 1), 'UniformOutput', false));
rt_end_array = cell2mat(cellfun(@(x) x(end), ric(:, 1), 'UniformOutput', false));
rt_tolerance = (rt_end_array - rt_start_array) ./ 5;
rt_intervals = [rt_start_array - rt_tolerance, rt_end_array + rt_tolerance];
[sorted_intervals, sort_idx] = sortrows(rt_intervals);
[categorized_intervals, categorized_indices] = categorize_intervals(sorted_intervals);
for idx_cat = 1:max(categorized_indices)
    % Extract the retention times and intensities of each IMP in the current category
    group_current_ric = ric(sort_idx(categorized_indices == idx_cat), :);
    group_current_imp_name = current_imp_name(sort_idx(categorized_indices == idx_cat));
    layout = CIMPQuantPlotter.getXicLegendLayoutConfig();
    plotXicGroupWithLayout(group_current_ric, total_xic, categorized_intervals(idx_cat, :), ...
        group_current_imp_name, fullfile(dir_save, [raw_name, '_', ...
        num2str(low_mz_bound), '-', num2str(high_mz_bound), '_+', ...
        num2str(selected_charge), '_', num2str(idx_cat)]), ...
        color_map, legend_map, layout);
end
end

function [categorized_intervals, categorized_indices] = categorize_intervals(intervals)
% Helper function to categorize retention time intervals
% input:
%   intervals (M x 2 double) minutes
%       retention time intervals
% output:
%   categorized_intervals (C x 2 double)
%       categorized retention time intervals; each row is [start_rt, end_rt] in minutes
%   categorized_indices (M x 1 double)
%       indices of the categories

categorized_intervals = zeros(0, 2); % Initialize an empty array to store categories
categorized_indices = zeros(size(intervals, 1), 1); % Initialize an array to store the category index of each interval
for i = 1:size(intervals, 1)
    % Extract the current time interval
    currentInterval = intervals(i, :);
    
    % Initialize a flag to indicate if the current interval has been categorized
    categorized = false;
    
    % Iterate over the existing categories
    for j = 1:size(categorized_intervals, 1)
        % Check if the current interval intersects with any interval in the category
        if is_intersecting(currentInterval, categorized_intervals(j, :))
            % If there is an intersection, add the current interval to this category
            categorized_intervals(j, 1) = min(currentInterval(1), categorized_intervals(j, 1));
            categorized_intervals(j, 2) = max(currentInterval(2), categorized_intervals(j, 2));
            categorized_indices(i) = j;
            categorized = true;
            break;
        end
    end
    
    % If the current interval does not intersect with any existing category, create a new category
    if ~categorized
        categorized_intervals(end+1, :) = currentInterval; %#ok<AGROW> 
        categorized_indices(i) = size(categorized_intervals, 1);
    end
end
end



function intersection = is_intersecting(interval1, interval2)
% Helper function to check if two time intervals intersect
% input:
%   interval1 (1 x 2 double) minutes
%   interval2 (1 x 2 double) minutes
% output:
%   intersection (1 x 1 logical)
intersection = (interval1(1) <= interval2(2)) && (interval1(2) >= interval2(1));
end
