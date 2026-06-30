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
    mz_tol = obj.m_spectrum.pre_mz * tolerance.value / 1e6;
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
set(gcf,'position',[50,50,900,600], 'color','white')
set(gca,'LooseInset',get(gca,'TightInset'));
set(gca,'Position',[0.12,0.11,0.84,0.82]);
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
if ~isfield(obj.m_peptides, 'label')
    return
end

raw_label = obj.m_peptides(peptide_index).label;
if isempty(raw_label)
    return
end

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
