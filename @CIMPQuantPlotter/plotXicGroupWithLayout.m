function diagnostics = plotXicGroupWithLayout(ric, total_xic, categorized_intervals, ...
    current_imp_name, file_base_path, color_map, legend_map, layout)
% Plot one XIC group with an adaptive right-side legend.
if nargin < 8 || isempty(layout)
    layout = CIMPQuantPlotter.getXicLegendLayoutConfig();
end

legend_x_fig = layout.left_margin_px + layout.axes_width_px + layout.axes_legend_gap_px;
if isfield(layout, 'figure_width_px') && ~isempty(layout.figure_width_px)
    fig_w = layout.figure_width_px;
else
    fig_w = legend_x_fig + layout.legend_min_width_px + layout.right_margin_px;
end
fig_h = layout.figure_height_px;
axes_height_px = fig_h - layout.bottom_margin_px - layout.top_margin_px;

f = figure('Visible', 'off', 'Units', 'pixels', 'Position', [50, 50, fig_w, fig_h], 'Color', 'white');
all_font_size = layout.all_font_size;
all_line_width = layout.all_line_width;
ax = axes('Parent', f, 'Units', 'pixels', 'FontSize', all_font_size);
hold(ax, 'on');

total_idx_start = find(total_xic{1} <= categorized_intervals(1), 1, 'last');
total_idx_end = find(total_xic{1} >= categorized_intervals(2), 1);
if isempty(total_idx_start) || total_idx_start < 1
    total_idx_start = 1;
end
if isempty(total_idx_end) || total_idx_end > numel(total_xic{1})
    total_idx_end = numel(total_xic{1});
end

plot(ax, total_xic{1}(total_idx_start:total_idx_end), total_xic{2}(total_idx_start:total_idx_end), ...
    'k', 'DisplayName', 'Total XIC', 'LineWidth', all_line_width);

plot_info = struct('x_data', {}, 'y_data', {}, 'legend_string', {}, 'raw_legend_string', {}, ...
    'color', {}, 'line_handle', {});
plot_count = 0;
string_to_color = @(str) hsv2rgb([mod(hash_string_positional(str), 360) / 360, 0.7, 0.9]);

for idx_imp = 1:size(ric, 1)
    if trapz(ric{idx_imp, 1}, ric{idx_imp, 2}) < 1e-6
        continue;
    end
    plot_count = plot_count + 1;
    plot_info(plot_count).x_data = ric{idx_imp, 1};
    plot_info(plot_count).y_data = ric{idx_imp, 2};
    is_in_legend_map = ~isempty(legend_map) && legend_map.isKey(current_imp_name{idx_imp});
    is_in_color_map = ~isempty(color_map) && color_map.isKey(current_imp_name{idx_imp});
    if is_in_legend_map
        plot_info(plot_count).legend_string = char(string(legend_map(current_imp_name{idx_imp})));
    else
        plot_info(plot_count).legend_string = char(string(current_imp_name{idx_imp}));
    end
    plot_info(plot_count).raw_legend_string = plot_info(plot_count).legend_string;
    if is_in_color_map
        plot_info(plot_count).color = color_map(current_imp_name{idx_imp});
    else
        plot_info(plot_count).color = string_to_color(current_imp_name{idx_imp});
    end
end

if plot_count > 0
    [~, sort_idx] = sort({plot_info.legend_string});
    plot_info = plot_info(sort_idx);
end

for idx_plot = 1:plot_count
    imp_xic_plot = plot(ax, plot_info(idx_plot).x_data, plot_info(idx_plot).y_data, ...
        'DisplayName', plot_info(idx_plot).legend_string, 'Color', plot_info(idx_plot).color, ...
        'LineWidth', all_line_width);
    plot_info(idx_plot).line_handle = imp_xic_plot;
end

xlabel(ax, 'Retention Time (min)', 'FontSize', all_font_size);
ylabel(ax, 'Intensity', 'FontSize', all_font_size);

ax_pos = [layout.left_margin_px, layout.bottom_margin_px, layout.axes_width_px, axes_height_px];
set(ax, 'Position', ax_pos);

line_chars = getXicLegendMaxLineChars(layout);
if isfield(layout, 'legend_min_line_chars') && ~isempty(layout.legend_min_line_chars)
    min_line_chars = layout.legend_min_line_chars;
else
    min_line_chars = 1;
end
min_export_margin_px = 8;
diagnostics = struct();
h_legend = [];
legend_col = legend_x_fig;
png_file = [file_base_path, '.png'];
final_legend_labels = {};
export_retry_count = 0;
max_export_retries = 3;

while true
    while true
        if ~isempty(h_legend) && isgraphics(h_legend)
            delete(h_legend);
        end

        if plot_count > 0
            prepared_labels = CIMPQuantPlotter.prepareXicLegendLabels( ...
                {plot_info.raw_legend_string}, layout, line_chars);
            for idx_plot = 1:plot_count
                plot_info(idx_plot).legend_string = prepared_labels{idx_plot};
                set(plot_info(idx_plot).line_handle, 'DisplayName', prepared_labels{idx_plot});
            end
            final_legend_labels = prepared_labels;
        end

        h_legend = legend(ax, 'show', 'Location', 'none');
        set(h_legend, 'FontSize', all_font_size, 'Box', 'off', 'Interpreter', 'none');
        drawnow;

        text_bounds = measureXicLegendTextBoundsFigurePx(ax, h_legend);
        print_width_scale = getLegendPrintWidthScale(layout);
        scaled_text_width = text_bounds.max_text_width * print_width_scale;
        max_allowed_text_width = layout.legend_max_width_px - layout.legend_padding_px - text_bounds.icon_width;
        if scaled_text_width > max_allowed_text_width && line_chars > min_line_chars
            line_chars = line_chars - 2;
            continue;
        end

        legend_width_px = min(max(text_bounds.icon_width + scaled_text_width + layout.legend_padding_px, ...
            layout.legend_min_width_px), layout.legend_max_width_px);
        content_height_px = max(text_bounds.content_height, 0);
        legend_height_px = content_height_px + layout.legend_padding_px;

        target_center_y = getAxesPlotCenterYFigurePx(ax, ax_pos);
        legend_y_fig = target_center_y - legend_height_px / 2;
        legend_y_fig = max(0, min(legend_y_fig, fig_h - legend_height_px));

        for idx_iter = 1:8
            % Legend Position with Units='pixels' is in figure coordinates, not axes-local.
            set(h_legend, 'Units', 'pixels', 'Position', ...
                [legend_x_fig, legend_y_fig, legend_width_px, legend_height_px]);
            drawnow;

            ink_bounds = measureXicLegendInkBoundsFigurePx(f, ax, layout, fig_w, fig_h);
            delta_y = target_center_y - ink_bounds.center_y;
            if abs(delta_y) < 1
                break;
            end
            legend_y_fig = legend_y_fig + delta_y;
            legend_y_fig = max(0, min(legend_y_fig, fig_h - legend_height_px));
        end

        set(f, 'Position', [50, 50, fig_w, fig_h]);
        set(ax, 'Position', ax_pos);
        drawnow;

        diagnostics = getXicLegendLayoutDiagnostics(f, ax, h_legend, layout);
        rendered_margin_px = measureRenderedLegendRightMarginPx(f, legend_col, fig_w);
        estimated_legend_right = legend_x_fig + legend_width_px;
        if diagnostics.is_legend_inside_figure && diagnostics.is_text_inside_figure && ...
                diagnostics.is_text_clear_of_axes && ...
                diagnostics.min_text_right_margin_px >= min_export_margin_px && ...
                estimated_legend_right <= fig_w - min_export_margin_px && ...
                rendered_margin_px >= min_export_margin_px
            break;
        end

        if line_chars <= min_line_chars
            close(f);
            error('CIMPQuantPlotter:LegendLayoutClipped', ...
                ['Legend layout diagnostics failed before export: %s. ', ...
                'Longest label: "%s". Legend budget: %.1f px. Figure width: %.1f px. ', ...
                'Increase figure_width_px/right legend budget or shorten legend labels.'], ...
                diagnostics.summary, getLongestLegendLabel(plot_info), layout.legend_max_width_px, fig_w);
        end
        line_chars = max(min_line_chars, line_chars - 4);
    end

    exportXicFigureWithLayout(f, file_base_path, layout, fig_w, fig_h);
    exported_margin_px = CIMPQuantPlotter.measureExportedLegendRightMarginPx(png_file, legend_col);
    if exported_margin_px >= min_export_margin_px
        assertXicExportIntegrity(file_base_path, layout, fig_w, fig_h);
        break;
    end

    export_retry_count = export_retry_count + 1;
    if line_chars <= min_line_chars || export_retry_count > max_export_retries
        close(f);
        error('CIMPQuantPlotter:LegendLayoutClipped', ...
            ['Exported PNG legend right margin too small: %d px. ', ...
            'Longest label: "%s". Legend budget: %.1f px. Figure width: %.1f px. ', ...
            'Increase figure_width_px/right legend budget or shorten legend labels.'], ...
            exported_margin_px, getLongestLegendLabel(plot_info), layout.legend_max_width_px, fig_w);
    end
    line_chars = max(min_line_chars, floor(line_chars / 2));
end

diagnostics.figure_position_px = [50, 50, fig_w, fig_h];
diagnostics.final_legend_labels = final_legend_labels;
diagnostics.final_legend_line_chars = line_chars;
diagnostics.rendered_legend_right_margin_px = measureRenderedLegendRightMarginPx(f, legend_col, fig_w);
close(f);
end

function scale = getLegendPrintWidthScale(layout)
if isfield(layout, 'legend_print_width_scale') && ~isempty(layout.legend_print_width_scale)
    scale = layout.legend_print_width_scale;
else
    scale = 1.0;
end
end

function hash_val = hash_string_positional(str)
if isempty(str)
    hash_val = 0;
    return;
end
char_codes = double(str);
hash1 = 0;
base = 31;
for i = 1:numel(char_codes)
    hash1 = mod(hash1 * base + char_codes(i), 2^32);
end
positions = 1:numel(char_codes);
hash2 = sum(char_codes .* positions .* 137);
signs = (-1) .^ (positions - 1);
hash3 = sum(char_codes .* signs .* positions);
hash_val = mod(hash1 + hash2 * 17 + abs(hash3) * 23, 360);
end

function label = getLongestLegendLabel(plot_info)
label = '';
max_len = 0;
for idx_plot = 1:numel(plot_info)
    current_label = char(string(plot_info(idx_plot).raw_legend_string));
    if numel(current_label) > max_len
        max_len = numel(current_label);
        label = current_label;
    end
end
end
