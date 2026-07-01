function plotMatch(obj, options)
% Plot the spectra and mark the matched peaks in the spectrum
% Input:
%   obj (CReviewPSM)
%       Review instance
%   options (optional struct)
%       max_labels_per_peak: maximum matched ion labels to show per peak
%       show_precursor: true to mark the precursor peak, false to show it
%           as a regular peak

if nargin < 2 || isempty(options)
    options = struct();
end
max_labels_per_peak = get_max_labels_per_peak(options);
show_precursor = get_show_precursor(options);
layout = get_msms_layout(options);
color_config = get_imp_color_config(obj, options);
label_y_offset = 1;

% Get relative intensity of all peaks
experimental_peaks = obj.m_spectrum.peaks;
experimental_peaks(:,2) = experimental_peaks(:,2)/max(experimental_peaks(:,2))*100;

% Separate the peak indexes into 3 group (unmatched, single matched and 
%   multiple matched), according to obj.m_all_match_ions:
%   [match_type, match_pos, charge, (*)expe_which, pep_which]
%   the fourth attribution
[matched_peaks, ~, matched_peaks_idxs] = unique(obj.m_all_match_ions(:,4),'rows');
count_idx = zeros(size(matched_peaks,1),1); % the number of appearance of each peak (row)
for idx = 1:length(matched_peaks_idxs)
    count_idx(matched_peaks_idxs(idx)) = count_idx(matched_peaks_idxs(idx))+1;
end
single_peak_indexes = matched_peaks(count_idx==1);
multi_peak_indexes = matched_peaks(count_idx>1);

% match precursor -> fragment ion
if obj.m_tolerance.is_ppm
    mz_tol = obj.m_spectrum.pre_mz * obj.m_tolerance.value / 1e6;
else
    mz_tol = obj.m_spectrum.pre_mz;
end
precursor_index = find(abs(experimental_peaks(:,1)-obj.m_spectrum.pre_mz)<mz_tol);
if length(precursor_index) > 1
    % The highest peak as precursor peak
    [~,precursor_index] = max(experimental_peaks(precursor_index,2));
end

% Plot 
figure
set(gcf,'position',[50,50,layout.figure_width_px,layout.figure_height_px], 'color','white')
ax = gca;
set(ax,'LooseInset',get(ax,'TightInset'));
set(ax,'Position',[0.12,0.11,layout.axes_width_fraction,0.82]);
hold on
ion_type = {'y','b'};

% Plot the precursor peak
if show_precursor
    X_precursor = experimental_peaks(precursor_index,1); % m/z of peaks
    Y_precursor = experimental_peaks(precursor_index,2); % intensity of peaks
    stem(X_precursor,Y_precursor,'|g','LineWidth',0.75); % stem plot
    text(X_precursor,Y_precursor+label_y_offset,['[M]',repmat('+',1,obj.m_spectrum.pre_charge)],'Rotation',90,'HorizontalAlignment','left')
else
    precursor_index = [];
end

% Plot the unmatched peaks
unmatched_peak_indexes = 1:size(experimental_peaks,1);
unmatched_peak_indexes([single_peak_indexes;multi_peak_indexes;precursor_index]) = [];
X_unmatched = experimental_peaks(unmatched_peak_indexes,1); % m/z of peaks
Y_unmatched = experimental_peaks(unmatched_peak_indexes,2); % intensity of peaks
stem(X_unmatched,Y_unmatched,'|','Color','#9E9E9E','LineWidth',0.75)            % stem plot


if color_config.is_enabled
    used_imp_labels = plot_colored_matched_peaks( ...
        obj, experimental_peaks, matched_peaks, color_config);

    for idx = 1:size(matched_peaks,1)
        selected_ions = obj.m_all_match_ions( ...
            find(obj.m_all_match_ions(:,4)==matched_peaks(idx)),:);
        ion_label = format_grouped_ion_labels(obj, selected_ions, ion_type, max_labels_per_peak);
        text(experimental_peaks(matched_peaks(idx),1), ...
            experimental_peaks(matched_peaks(idx),2)+label_y_offset, ...
            ion_label,'Rotation',90,'HorizontalAlignment','left')
    end

    if color_config.show_legend
        draw_imp_color_legend(used_imp_labels, color_config);
    end
else
    % Plot the single matched peaks

    X_single = experimental_peaks(single_peak_indexes,1);
    Y_single = experimental_peaks(single_peak_indexes,2);
    stem(X_single,Y_single,'|r','LineWidth',1.5)
    for idx = 1:size(single_peak_indexes,1)  % labels of single matched peaks
        selected_ions = obj.m_all_match_ions(...
            find(obj.m_all_match_ions(:,4)==single_peak_indexes(idx)),:);
        ion_label = format_grouped_ion_labels(obj, selected_ions, ion_type, max_labels_per_peak);
        text(X_single(idx),Y_single(idx)+label_y_offset,ion_label,'Rotation',90,'HorizontalAlignment','left')
    %     text(X_single(idx),Y_single(idx)+0.01,ion_label,'HorizontalAlignment','center')
    end

    % Plot the multi matched peaks
    hold on
    X_single = experimental_peaks(multi_peak_indexes,1);
    Y_single = experimental_peaks(multi_peak_indexes,2);
    stem(X_single,Y_single,'|b','LineWidth',1.5)
    for idx = 1:size(multi_peak_indexes,1)  % labels of multi matched peaks
        selected_ions = obj.m_all_match_ions(...
            find(obj.m_all_match_ions(:,4)==multi_peak_indexes(idx)),:);
        ion_label = format_grouped_ion_labels(obj, selected_ions, ion_type, max_labels_per_peak);
        text(X_single(idx),Y_single(idx)+label_y_offset,ion_label,'Rotation',90,'HorizontalAlignment','left')
    end
end
xlabel('m/z')
ylim([0, max(115, max(experimental_peaks(:,2)) + 15)]);
set(gca, 'YTick', 0:20:100, ...
    'YTickLabel', {'0', '20', '40', '60', '80', '100'});
draw_custom_y_axis(gca, 0:20:100);
% saveas(gcf,'test.svg') % save into vector graph
end

function max_labels_per_peak = get_max_labels_per_peak(options)
% Return the configured per-peak label limit, preserving the legacy default.

max_labels_per_peak = 2;
if isfield(options, 'max_labels_per_peak') && ~isempty(options.max_labels_per_peak)
    max_labels_per_peak = options.max_labels_per_peak;
end
if ~isscalar(max_labels_per_peak) || max_labels_per_peak < 1 || ...
        max_labels_per_peak ~= floor(max_labels_per_peak)
    error('CReviewPSM:InvalidMaxLabelsPerPeak', ...
        'options.max_labels_per_peak must be a positive integer scalar.');
end
end

function show_precursor = get_show_precursor(options)
% Return whether precursor peaks should receive a special marker.

show_precursor = true;
if isfield(options, 'show_precursor') && ~isempty(options.show_precursor)
    show_precursor = options.show_precursor;
end
if ~islogical(show_precursor) || ~isscalar(show_precursor)
    error('CReviewPSM:InvalidShowPrecursor', ...
        'options.show_precursor must be a logical scalar.');
end
end

function layout = get_msms_layout(options)
% Return figure and axes geometry for the MS/MS matching plot.

layout = struct( ...
    'figure_width_px', 900, ...
    'figure_height_px', 600, ...
    'axes_width_fraction', 0.84);
if ~isfield(options, 'msms_layout') || isempty(options.msms_layout)
    return
end
if ~isstruct(options.msms_layout) || ~isscalar(options.msms_layout)
    error('CReviewPSM:InvalidMsmsLayout', ...
        'options.msms_layout must be a scalar struct.');
end
requested_layout = options.msms_layout;
layout = copy_positive_scalar_field(layout, requested_layout, 'figure_width_px');
layout = copy_positive_scalar_field(layout, requested_layout, 'figure_height_px');
layout = copy_axes_width_fraction(layout, requested_layout);
end

function layout = copy_positive_scalar_field(layout, requested_layout, field_name)
if ~isfield(requested_layout, field_name) || isempty(requested_layout.(field_name))
    return
end
value = requested_layout.(field_name);
if ~isnumeric(value) || ~isscalar(value) || ~isfinite(value) || value <= 0
    error('CReviewPSM:InvalidMsmsLayout', ...
        'options.msms_layout.%s must be a positive numeric scalar.', field_name);
end
layout.(field_name) = value;
end

function layout = copy_axes_width_fraction(layout, requested_layout)
field_name = 'axes_width_fraction';
if ~isfield(requested_layout, field_name) || isempty(requested_layout.(field_name))
    return
end
value = requested_layout.(field_name);
if ~isnumeric(value) || ~isscalar(value) || ~isfinite(value) || ...
        value <= 0 || value > 0.84
    error('CReviewPSM:InvalidMsmsLayout', ...
        'options.msms_layout.axes_width_fraction must be in the range (0, 0.84].');
end
layout.axes_width_fraction = value;
end

function color_config = get_imp_color_config(obj, options)
% Parse optional IMP color/proportion maps and decide whether coloring can run.

color_config = struct( ...
    'is_enabled', false, ...
    'colors', [], ...
    'proportions', [], ...
    'show_legend', get_show_imp_color_legend(options));

has_colors = isfield(options, 'imp_colors') && ~isempty(options.imp_colors);
has_proportions = isfield(options, 'imp_proportions') && ~isempty(options.imp_proportions);
if has_colors
    validate_color_map(options.imp_colors);
end
if has_proportions
    validate_proportion_map(options.imp_proportions);
end
if ~(has_colors && has_proportions)
    return
end

matched_peptide_indexes = unique(obj.m_all_match_ions(:,5));
if isempty(matched_peptide_indexes)
    return
end
for idx_peptide = 1:numel(matched_peptide_indexes)
    peptide_index = matched_peptide_indexes(idx_peptide);
    if ~has_explicit_peptide_label(obj, peptide_index)
        return
    end
    peptide_label = get_peptide_label(obj, peptide_index);
    if ~options.imp_colors.isKey(peptide_label) || ...
            ~options.imp_proportions.isKey(peptide_label)
        return
    end
end
validate_matched_imp_proportion_sum(obj, matched_peptide_indexes, options.imp_proportions);
validate_no_zero_proportion_matched_peaks(obj, options.imp_proportions);

color_config.is_enabled = true;
color_config.colors = options.imp_colors;
color_config.proportions = options.imp_proportions;
end

function validate_matched_imp_proportion_sum(obj, matched_peptide_indexes, proportion_map)
matched_sum = 0;
for idx_peptide = 1:numel(matched_peptide_indexes)
    peptide_label = get_peptide_label(obj, matched_peptide_indexes(idx_peptide));
    matched_sum = matched_sum + proportion_map(peptide_label);
end
if matched_sum <= 0
    error('CReviewPSM:InvalidImpProportion', ...
        'The matched IMP proportions for a colored plot must have a positive sum.');
end
end

function validate_no_zero_proportion_matched_peaks(obj, proportion_map)
matched_peak_indexes = unique(obj.m_all_match_ions(:, 4));
for idx_peak = 1:numel(matched_peak_indexes)
    peak_index = matched_peak_indexes(idx_peak);
    selected_ions = obj.m_all_match_ions(obj.m_all_match_ions(:,4)==peak_index, :);
    peak_imp_labels = get_unique_ion_imp_labels(obj, selected_ions);
    peak_proportion_sum = 0;
    for idx_label = 1:numel(peak_imp_labels)
        peak_proportion_sum = peak_proportion_sum + proportion_map(peak_imp_labels{idx_label});
    end
    if peak_proportion_sum <= 0
        error('CReviewPSM:ZeroImpProportionMatchedPeak', ...
            'Matched peak %d only maps to zero-proportion IMP labels: %s.', ...
            peak_index, strjoin(peak_imp_labels, ','));
    end
end
end

function show_legend = get_show_imp_color_legend(options)
% Return whether the optional IMP color legend should be drawn.

show_legend = false;
if isfield(options, 'show_imp_color_legend') && ~isempty(options.show_imp_color_legend)
    show_legend = options.show_imp_color_legend;
end
if ~islogical(show_legend) || ~isscalar(show_legend)
    error('CReviewPSM:InvalidShowImpColorLegend', ...
        'options.show_imp_color_legend must be a logical scalar.');
end
end

function validate_color_map(color_map)
if ~isa(color_map, 'containers.Map')
    error('CReviewPSM:InvalidImpColors', ...
        'options.imp_colors must be a containers.Map.');
end
color_keys = keys(color_map);
for idx_key = 1:numel(color_keys)
    color_value = color_map(color_keys{idx_key});
    if ~isnumeric(color_value) || ~isequal(size(color_value), [1, 3]) || ...
            any(~isfinite(color_value)) || any(color_value < 0) || any(color_value > 1)
        error('CReviewPSM:InvalidImpColor', ...
            'Each options.imp_colors value must be a finite 1x3 RGB vector in [0, 1].');
    end
end
end

function validate_proportion_map(proportion_map)
if ~isa(proportion_map, 'containers.Map')
    error('CReviewPSM:InvalidImpProportions', ...
        'options.imp_proportions must be a containers.Map.');
end
proportion_keys = keys(proportion_map);
for idx_key = 1:numel(proportion_keys)
    proportion_value = proportion_map(proportion_keys{idx_key});
    if ~isnumeric(proportion_value) || ~isscalar(proportion_value) || ...
            ~isfinite(proportion_value) || proportion_value < 0
        error('CReviewPSM:InvalidImpProportion', ...
            'Each options.imp_proportions value must be a finite non-negative scalar.');
    end
end
end

function used_imp_labels = plot_colored_matched_peaks(obj, experimental_peaks, ...
    matched_peak_indexes, color_config)
% Draw matched peaks using IMP colors; shared peaks are height-stacked.

used_imp_labels = {};
for idx_peak = 1:size(matched_peak_indexes, 1)
    peak_index = matched_peak_indexes(idx_peak);
    selected_ions = obj.m_all_match_ions(obj.m_all_match_ions(:,4)==peak_index, :);
    peak_imp_labels = get_unique_ion_imp_labels(obj, selected_ions);
    peak_proportions = zeros(1, numel(peak_imp_labels));
    for idx_label = 1:numel(peak_imp_labels)
        peak_proportions(idx_label) = color_config.proportions(peak_imp_labels{idx_label});
    end
    proportion_sum = sum(peak_proportions);
    if proportion_sum <= 0
        error('CReviewPSM:ZeroImpProportionMatchedPeak', ...
            'Matched peak %d only maps to zero-proportion IMP labels: %s.', ...
            peak_index, strjoin(peak_imp_labels, ','));
    end
    peak_proportions = peak_proportions ./ proportion_sum;

    x_value = experimental_peaks(peak_index, 1);
    y_value = experimental_peaks(peak_index, 2);
    y_start = 0;
    for idx_label = 1:numel(peak_imp_labels)
        y_end = y_start + y_value * peak_proportions(idx_label);
        color_value = color_config.colors(peak_imp_labels{idx_label});
        if numel(peak_imp_labels) == 1
            object_tag = 'CReviewPSMImpColoredStem';
        else
            object_tag = 'CReviewPSMImpColoredSegment';
        end
        line([x_value, x_value], [y_start, y_end], ...
            'Color', color_value, ...
            'LineWidth', 1.5, ...
            'Tag', object_tag);
        y_start = y_end;
    end
    used_imp_labels = append_new_labels(used_imp_labels, peak_imp_labels);
end
end

function peak_imp_labels = get_unique_ion_imp_labels(obj, selected_ions)
% Return unique peptide labels for all ions matched to one experimental peak.

peak_imp_labels = {};
for idx_ion = 1:size(selected_ions, 1)
    peptide_label = get_peptide_label(obj, selected_ions(idx_ion, 5));
    if ~any(strcmp(peak_imp_labels, peptide_label))
        peak_imp_labels{end + 1} = peptide_label; %#ok<AGROW>
    end
end
end

function combined_labels = append_new_labels(combined_labels, labels_to_add)
for idx_label = 1:numel(labels_to_add)
    if ~any(strcmp(combined_labels, labels_to_add{idx_label}))
        combined_labels{end + 1} = labels_to_add{idx_label}; %#ok<AGROW>
    end
end
end

function draw_imp_color_legend(used_imp_labels, color_config)
% Draw an outside legend for IMP color identities.

if isempty(used_imp_labels)
    return
end
legend_handles = gobjects(1, numel(used_imp_labels));
for idx_label = 1:numel(used_imp_labels)
    legend_handles(idx_label) = plot(nan, nan, '-', ...
        'Color', color_config.colors(used_imp_labels{idx_label}), ...
        'LineWidth', 1.5, ...
        'DisplayName', used_imp_labels{idx_label}, ...
        'Tag', 'CReviewPSMImpColorLegendHandle');
end
legend_handle = legend(legend_handles, used_imp_labels, ...
    'Location', 'eastoutside', ...
    'Interpreter', 'none');
legend_handle.Tag = 'CReviewPSMImpColorLegend';
end

function draw_custom_y_axis(ax, tick_values)
% Draw a y axis that visually stops at 100 while ylim leaves label headroom.

x_limits = xlim(ax);
x_left = x_limits(1);
x_span = diff(x_limits);
axis_color = [0.15, 0.15, 0.15];
tick_length = 0.010 * x_span;
tick_label_gap = 0.018 * x_span;
axis_label_gap = 0.055 * x_span;
font_size = get(ax, 'FontSize');
font_name = get(ax, 'FontName');

ax.YAxis.Visible = 'off';
set(ax, 'YTick', [], 'YTickLabel', {});

line(ax, [x_left, x_left], [0, 100], ...
    'Color', axis_color, ...
    'LineWidth', 0.6667, ...
    'Clipping', 'off', ...
    'HandleVisibility', 'off', ...
    'Tag', 'CReviewPSMCustomYAxisLine');

for idx_tick = 1:numel(tick_values)
    tick_value = tick_values(idx_tick);
    line(ax, [x_left, x_left + tick_length], [tick_value, tick_value], ...
        'Color', axis_color, ...
        'LineWidth', 0.6667, ...
        'Clipping', 'off', ...
        'HandleVisibility', 'off', ...
        'Tag', 'CReviewPSMCustomYAxisTick');
    text(ax, x_left - tick_label_gap, tick_value, num2str(tick_value), ...
        'HorizontalAlignment', 'right', ...
        'VerticalAlignment', 'middle', ...
        'FontSize', font_size, ...
        'FontName', font_name, ...
        'Color', axis_color, ...
        'Clipping', 'off', ...
        'Tag', 'CReviewPSMCustomYAxisTickLabel');
end

text(ax, x_left - axis_label_gap, 50, 'Relative intensity (%)', ...
    'Rotation', 90, ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', ...
    'FontSize', font_size + 1, ...
    'FontName', font_name, ...
    'Color', axis_color, ...
    'Clipping', 'off', ...
    'Tag', 'CReviewPSMCustomYAxisLabel');
end

function ion_label = format_grouped_ion_labels(obj, selected_ions, ion_type, max_labels_per_peak)
% Format matched ions by grouping peptide labels that share the same ion.

group_labels = {};
group_keys = zeros(0, 3);
for idx_ion = 1:size(selected_ions, 1)
    selected_ion = selected_ions(idx_ion, :);
    ion_key = selected_ion(1:3);
    group_idx = find(ismember(group_keys, ion_key, 'rows'), 1);
    if isempty(group_idx)
        group_keys(end + 1, :) = ion_key; %#ok<AGROW>
        group_labels{end + 1} = {get_peptide_label(obj, selected_ion(5))}; %#ok<AGROW>
    else
        peptide_label = get_peptide_label(obj, selected_ion(5));
        if ~any(strcmp(group_labels{group_idx}, peptide_label))
            group_labels{group_idx}{end + 1} = peptide_label;
        end
    end
end

label_count = min(size(group_keys, 1), max_labels_per_peak);
label_parts = cell(1, label_count);
for idx_group = 1:label_count
    ion_text = format_ion_text(group_keys(idx_group, :), ion_type);
    label_parts{idx_group} = [strjoin(group_labels{idx_group}, ','), ':', ion_text];
end
ion_label = strjoin(label_parts, ';');
if size(group_keys, 1) > max_labels_per_peak
    ion_label = [ion_label, ';...'];
end
end

function ion_text = format_ion_text(ion_key, ion_type)
% Format one ion identity without peptide ownership.

ion_text = [ion_type{ion_key(1)}, num2str(ion_key(2)), ...
    repmat('+',1,ion_key(3))];
end

function peptide_label = get_peptide_label(obj, peptide_index)
% Return the optional peptide display label, or fall back to the local index.

peptide_label = num2str(peptide_index);
if ~has_explicit_peptide_label(obj, peptide_index)
    return
end

raw_label = obj.m_peptides(peptide_index).label;
if ischar(raw_label)
    peptide_label = raw_label;
elseif isstring(raw_label) && isscalar(raw_label)
    peptide_label = char(raw_label);
elseif isnumeric(raw_label) && isscalar(raw_label)
    peptide_label = num2str(raw_label);
else
    peptide_label = char(string(raw_label));
end
end

function has_label = has_explicit_peptide_label(obj, peptide_index)
% Return true only when a peptide has a non-empty explicit display label.

has_label = false;
if ~isfield(obj.m_peptides, 'label')
    return
end
raw_label = obj.m_peptides(peptide_index).label;
if isempty(raw_label)
    return
end
has_label = true;
end
