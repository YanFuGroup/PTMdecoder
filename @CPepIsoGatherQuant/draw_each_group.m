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

% System error in saving retention time
eps_rt_print = 1e-6;

num_iso = size(current_ratioMatrix,2);
% auxic = zeros(num_iso, 1);
ric = cell(num_iso, 2); % IMP -> rt_grid, intensity, total intensity
rt_error_tol = 1; % RT match tolerance, choose 1 arbitrarily
% A vector showing is needed to skip this IMP.
%   Cannot delete because other filter can also change this vector
is_skip_vec = cellfun(@isempty,current_iso_rt_range);

% Sort MS1 signal (pair of retention time and intensity) by time
sort_rts = [(1:length(current_rts))',current_rts];
sort_rts = sortrows(sort_rts,2);% Sort in ascending order by retention time column
sort_idx = sort_rts(:,1);
sort_rts = sort_rts(:,2);
sort_inten = current_inten(sort_idx);
sort_ratioMatrix = current_ratioMatrix(sort_idx,:);% Rearrange the matrix in chronological order

% Sort and denoise using a relative abundance threshold method
sort_inten = smooth(sort_inten,0.05,'loess');%0.05,'loess'
maxInten = max(sort_inten);
tmp = sort_inten<0.05*maxInten;% Find results where intensity is less than 0.05 of the maximum value and discard them
sort_inten(tmp) = []; %#ok<NASGU>
sort_rts(tmp) = [];
sort_ratioMatrix(tmp,:) = [];

if (obj.hasMinRows(sort_ratioMatrix, obj.m_minMSMSnum) == false)
    % If the ratio matrix has less than 3 rows, skip this group
    return;
end

% find the XIC filtered by m/z of base peak [low_mz_bound, high_mz_bound]
%   and -1,+0,+1,+2,+3
% MS1_index (scan, retention time, peak number, baseline, injection time)
% MS1_peaks (m/z, intensity)
MS1_index = obj.m_cMs12DatasetIO.m_mapNameMS1Index(erase(raw_name,'.mgf'));
MS1_peaks = obj.m_cMs12DatasetIO.m_mapNameMS1Peaks(erase(raw_name,'.mgf'));
rt_grid = MS1_index(:,2);
isotope_num = [-1,0,1,2,3,4];
intensity = zeros(size(MS1_index,1),length(isotope_num)); % retention time -> intensity, XIC
% record the isotopic XIC
for idx_iso = 1:length(isotope_num)
    idxs_target_peaks = find(MS1_peaks(:,1)>low_mz_bound+isotope_num(idx_iso)*CConstant.unitdiff/selected_charge...
        & MS1_peaks(:,1)<high_mz_bound+isotope_num(idx_iso)*CConstant.unitdiff/selected_charge);
    for idx_itp = 1:length(idxs_target_peaks)
        intensity(find(MS1_index(:,3)>idxs_target_peaks(idx_itp),1),idx_iso) = ...
            MS1_peaks(idxs_target_peaks(idx_itp),2);
    end
end
% filter with two criteria:
% 1. the intensity of -1 peak should not greater than monoisotopic
% 2. the intensity of isotopic cluster peak should be enough similar with
%   the IPV matrix.
for idx_inten = 1:size(intensity,1)
    if intensity(idx_inten,1)>intensity(idx_inten,2) % the first criterion
        intensity(idx_inten,:) = 0;
    elseif ~any(intensity(idx_inten,:))
        continue;
    elseif 1-pdist([CConstant.IPV(int64((high_mz_bound+low_mz_bound)/2),:);...
            intensity(idx_inten,2:end)],'cosine') < 0.6
        intensity(idx_inten,:) = 0;
    end
end
intensity = intensity(:,2);
smoothed_intensity = smoothdata(intensity,'movmean',5);
% smoothed_intensity = smoothdata(intensity,'gaussian',1,'SamplePoints',rt_grid);

% Extract the rt bound of XIC peak
final_XIC_peak_for_IMP = repmat(struct('left_bound',0,'right_bound',0), num_iso, 1);
max_label = zeros(num_iso,1);
for idx_iso = 1:num_iso
    % Check if need to consider this IMP
    if is_skip_vec(idx_iso)
        continue;
    end

    % Record all of the check labels
    check_labels = zeros(length(current_iso_rt_range{idx_iso}),1);
    for idx_peak = 1:length(current_iso_rt_range{idx_iso})
        check_labels(idx_peak) = current_iso_rt_range{idx_iso}(idx_peak).check_label;
    end
    % Find the peak with max check label (the first of max peaks)
    [max_label(idx_iso), idx_max] = max(check_labels);
    if max_label(idx_iso) == 0
        % If all check labels are zero, skip this IMP later on
        is_skip_vec(idx_iso)=true;
        continue;
    end
    final_XIC_peak_for_IMP(idx_iso).left_bound = current_iso_rt_range{idx_iso}(idx_max).rt_start;
    final_XIC_peak_for_IMP(idx_iso).right_bound = current_iso_rt_range{idx_iso}(idx_max).rt_end;
end

% Calculate the ratio on each XIC points using kernel method, and normalize
% using Nadaraya-Waston kernel averaging method
esti_ratio = zeros(length(rt_grid),num_iso);
% % bandwidth = (4/(3*size(sort_rts,1)))^0.2*std(sort_rts);
for idx_iso = 1:num_iso
    % Check if need to consider this IMP
    if is_skip_vec(idx_iso)
        continue;
    end

    % Gaussian kernel function
    Ker_Gaussian = @(u) (1/sqrt(2*pi))*exp(-0.5*u.^2);
    % Epanechnikov kernel function
    %     Ker_Epanechnikov = @(u) (3/4)*(1-u.^2).*(abs(u)<=1);
    
    % Retention time bound of current IMP
    cur_iso_rt_left = final_XIC_peak_for_IMP(idx_iso).left_bound;
    cur_iso_rt_right = final_XIC_peak_for_IMP(idx_iso).right_bound;
    [rt_diff, idx_rt_left] = min(abs(rt_grid-cur_iso_rt_left));
    if rt_diff > rt_error_tol
        error(['Cannot find the spectra on the specified retention time: ', ...
            num2str(cur_iso_rt_left)]);
    end
    [rt_diff, idx_rt_right] = min(abs(rt_grid-cur_iso_rt_right));
    if rt_diff > rt_error_tol
        error(['Cannot find the spectra on the specified retention time: ', ...
            num2str(cur_iso_rt_right)]);
    end

    % Collect all of the rts states within current XIC peak
    idxs_ident_rt = sort_rts>=cur_iso_rt_left-eps_rt_print & sort_rts<=cur_iso_rt_right+eps_rt_print;
    rts_current = sort_rts(idxs_ident_rt);

    % Calculate bandwidth and weights for each XIC peak
    bandwidth = (4/(3*size(rts_current,1)))^0.2*std(rts_current);
    weights = zeros(idx_rt_right-idx_rt_left+1, length(rts_current));
    for idx_PSM = 1:length(rts_current)
        if bandwidth == 0
            break;
        end
        % set kernel weights
        weights(:,idx_PSM) = Ker_Gaussian(...
            (rt_grid(idx_rt_left:idx_rt_right)-rts_current(idx_PSM))/bandwidth);
    end
    % Check if there are nearly no weights in some retention time for all
    %   IMP, or the bandwidth is just zero
    if bandwidth == 0 || any(all(weights<1e-15/length(sort_rts),2))
        bandwidth = min(cur_iso_rt_right-cur_iso_rt_left,1);
        for idx_PSM = 1:length(rts_current)
            weights(:,idx_PSM) = Ker_Gaussian(...
                (rt_grid(idx_rt_left:idx_rt_right)-rts_current(idx_PSM))/bandwidth);
        end
    end
    % Calculate the ratio using normalized weights and ratioMatrix
    esti_ratio(idx_rt_left:idx_rt_right,idx_iso) = ...
        (weights * sort_ratioMatrix(idxs_ident_rt,idx_iso))./(sum(weights,2)+eps);
end
% normalize the ratio in every available retention time
esti_ratio = esti_ratio./(sum(esti_ratio,2)+eps);

% Hypothesis: The retention time bounds revised manually cannot be wrong.
%   No need to recognize and delete the low abundance IMP in revised file.
% Filter low abundance IMP according to the relative AUXIC
% intensityMatrix = esti_ratio.*smoothed_intensity;
% for i_Xp = 1:length(XIC_peaks)
%     area_filter = zeros(num_iso,1);
%     for idx_iso = 1:num_iso
%         area_filter(idx_iso) = trapz(rt_grid(XIC_peaks(i_Xp).left_bound:XIC_peaks(i_Xp).right_bound),...
%             intensityMatrix(XIC_peaks(i_Xp).left_bound:XIC_peaks(i_Xp).right_bound,idx_iso));
%     end
%     % If the AUC of the peak of an IMP is less than 10% of the maximum,
%     % filter, remove the proportion.
%     esti_ratio(XIC_peaks(i_Xp).left_bound:XIC_peaks(i_Xp).right_bound,area_filter<max(area_filter)*obj.m_resFilterThres) = 0;
% end
% esti_ratio = esti_ratio./(sum(esti_ratio,2)+eps);

% Requantification using revised RT
% plot(rt_grid, smoothed_intensity, 'k', 'DisplayName','Total XIC');
total_xic = {rt_grid, smoothed_intensity};
intensityMatrix = esti_ratio.*smoothed_intensity;
rt_bound = repmat(struct('start',0,'end',0), num_iso, 1);
for idx_iso = 1:num_iso
    % Check if need to consider this IMP
    if is_skip_vec(idx_iso)
        continue;
    end

    % Get the final rt bound
    rt_bound(idx_iso).start = final_XIC_peak_for_IMP(idx_iso).left_bound;
    rt_bound(idx_iso).end = final_XIC_peak_for_IMP(idx_iso).right_bound;
    [rt_diff, final_rt_start] = min(abs(rt_grid-rt_bound(idx_iso).start));
    if rt_diff > rt_error_tol
        error(['Cannot find the spectra on the specified retention time: ', ...
            num2str(rt_bound(idx_iso).start)]);
    end
    if final_rt_start ~= 1
        final_rt_start = final_rt_start - 1;
    end
    [rt_diff, final_rt_end] = min(abs(rt_grid-rt_bound(idx_iso).end));
    if rt_diff > rt_error_tol
        error(['Cannot find the spectra on the specified retention time: ', ...
            num2str(rt_bound(idx_iso).end)]);
    end
    if final_rt_end ~= length(rt_grid)
        final_rt_end = final_rt_end + 1;
    end
    % auxic(idx_iso,1) = trapz(rt_grid(final_rt_start:final_rt_end),...
    %     [0;intensityMatrix(final_rt_start+1:final_rt_end-1,idx_iso);0])*60;
    % total_temp = trapz(rt_grid(final_rt_start:final_rt_end),...
    %     [0;smoothed_intensity(final_rt_start+1:final_rt_end-1);0])*60;
    % if auxic(idx_iso,1) < total_temp*0.01
    %     continue;
    % end
    ric{idx_iso,1} = rt_grid(final_rt_start:final_rt_end);
    ric{idx_iso,2} = [0;intensityMatrix(final_rt_start+1:final_rt_end-1,idx_iso);0];
    % ric{idx_iso,3} = [0;smoothed_intensity(final_rt_start+1:final_rt_end-1);0];
    % plot(rt_grid(final_rt_start:final_rt_end),[0;intensityMatrix(final_rt_start+1:final_rt_end-1,idx_iso);0],...
    %     'DisplayName', ['RIC of ',current_iso_name{idx_iso}]);
end

% Check if the output directory exists
if ~exist(dir_save,'dir')
    mkdir(dir_save);
end

% Plot the total XIC and RIC
plot_xics(ric, current_iso_name, total_xic, dir_save, raw_name, low_mz_bound, high_mz_bound, selected_charge, color_map, legend_map);
end



function plot_xics(ric, current_iso_name, total_xic, dir_save, raw_name, low_mz_bound, high_mz_bound, selected_charge, color_map, legend_map)
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
f = figure('visible','off');
set(gcf,'position',[50, 50, 600, 200], 'color','white')
% set(gcf,'position',[50, 50, 300, 200], 'color','white')
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
for idx_iso = 1:size(ric, 1)
    if trapz(ric{idx_iso, 1}, ric{idx_iso, 2}) < 1e-6
        continue;
    end
    
    plot_count = plot_count + 1;
    plot_info(plot_count).x_data = ric{idx_iso, 1};
    plot_info(plot_count).y_data = ric{idx_iso, 2};
    
    is_in_legend_map = ~isempty(legend_map) && legend_map.isKey(current_iso_name{idx_iso});
    is_in_color_map = ~isempty(color_map) && color_map.isKey(current_iso_name{idx_iso});

    if is_in_legend_map
        plot_info(plot_count).legend_string = ['XIC of ',legend_map(current_iso_name{idx_iso})];
    else
        legend_string = ['XIC of ',current_iso_name{idx_iso}];
        legend_string = strrep(legend_string, '_', '\_');
        legend_string = strrep(legend_string, '{', '\{');
        legend_string = strrep(legend_string, '}', '\}');
        plot_info(plot_count).legend_string = legend_string;
    end

    if is_in_color_map
        plot_info(plot_count).color = color_map(current_iso_name{idx_iso});
    else
        % Generate hash-based color from string
        plot_info(plot_count).color = string_to_color(current_iso_name{idx_iso});
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

h_legend = legend('show', 'Location', 'northwest');
set(h_legend, 'FontSize', all_font_size);
saveas(f, file_name);
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