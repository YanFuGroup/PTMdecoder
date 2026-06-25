function bounds = measureXicLegendRenderedBoundsFigurePx(ax, h_legend)
% Measure rendered legend content bounds in figure pixel coordinates.
% Uses a top-down proxy model because MATLAB legend objects do not expose
% Text children in batch/off-screen rendering.
set(ax, 'Units', 'pixels');
set(h_legend, 'Units', 'pixels');
legend_pos_fig = get(h_legend, 'Position');

bounds = measure_top_down_proxy(ax, legend_pos_fig, legend_pos_fig, h_legend);
end

function bounds = measure_top_down_proxy(ax, legend_pos_fig, legend_pos_ax, h_legend)
strings = h_legend.String;
if ischar(strings)
    strings = cellstr(strings);
elseif isstring(strings)
    strings = cellstr(strings);
end

icon_width = 30;
if isprop(h_legend, 'ItemTokenSize')
    icon_width = h_legend.ItemTokenSize(1);
end
icon_height = icon_width;
if isprop(h_legend, 'ItemTokenSize') && numel(h_legend.ItemTokenSize) >= 2
    icon_height = h_legend.ItemTokenSize(2);
end
font_size = h_legend.FontSize;
font_name = h_legend.FontName;
entry_gap_px = 6;

legend_top = legend_pos_fig(2) + legend_pos_fig(4);
cursor_top = legend_top;
text_extents = zeros(0, 4);
line_extents = zeros(0, 4);

for idx_string = 1:numel(strings)
    label_lines = splitXicLegendLabelLines(strings{idx_string});
    entry_line_heights = zeros(numel(label_lines), 1);
    for idx_line = 1:numel(label_lines)
        text_extent = measureXicLegendLabelTextExtentPx( ...
            ax, label_lines{idx_line}, font_size, font_name);
        entry_line_heights(idx_line) = text_extent(4);
    end
    entry_height = max(sum(entry_line_heights), icon_height);
    entry_bottom = cursor_top - entry_height;
    icon_center_y = entry_bottom + entry_height / 2;
    line_extents(end + 1, :) = [legend_pos_fig(1), icon_center_y - icon_height / 2, icon_width, icon_height]; %#ok<AGROW>

    line_cursor_top = cursor_top;
    for idx_line = 1:numel(label_lines)
        line_height = entry_line_heights(idx_line);
        line_bottom = line_cursor_top - line_height;
        text_extent = measureXicLegendLabelTextExtentPx( ...
            ax, label_lines{idx_line}, font_size, font_name);
        text_extents(end + 1, :) = [legend_pos_fig(1) + icon_width, line_bottom, text_extent(3), line_height]; %#ok<AGROW>
        line_cursor_top = line_bottom;
    end

    cursor_top = entry_bottom;
    if idx_string < numel(strings)
        cursor_top = cursor_top - entry_gap_px;
    end
end

all_extents = [text_extents; line_extents];
if isempty(all_extents)
    bounds = empty_bounds(legend_pos_fig);
    return;
end

min_left = min(all_extents(:, 1));
max_right = max(all_extents(:, 1) + all_extents(:, 3));
min_bottom = min(all_extents(:, 2));
max_top = max(all_extents(:, 2) + all_extents(:, 4));
max_text_width = 0;
if ~isempty(text_extents)
    max_text_width = max(text_extents(:, 3));
end

bounds = pack_bounds(min_left, max_right, min_bottom, max_top, text_extents, ...
    legend_pos_fig, icon_width);
bounds.max_text_width = max_text_width;
end

function bounds = pack_bounds(min_left, max_right, min_bottom, max_top, text_extents, ...
    legend_pos_fig, icon_width)
content_height = max_top - min_bottom;
if isempty(text_extents)
    max_text_width = 0;
else
    max_text_width = max(text_extents(:, 3));
end
bounds = struct('min_left', min_left, 'max_right', max_right, ...
    'min_bottom', min_bottom, 'max_top', max_top, ...
    'content_height', content_height, ...
    'text_extents', text_extents, 'legend_position_px', legend_pos_fig, ...
    'icon_width', icon_width, 'max_text_width', max_text_width);
end

function bounds = empty_bounds(legend_pos_fig)
bounds = struct('min_left', legend_pos_fig(1), 'max_right', legend_pos_fig(1) + legend_pos_fig(3), ...
    'min_bottom', legend_pos_fig(2), 'max_top', legend_pos_fig(2) + legend_pos_fig(4), ...
    'content_height', legend_pos_fig(4), ...
    'text_extents', zeros(0, 4), 'legend_position_px', legend_pos_fig, ...
    'icon_width', 0, 'max_text_width', 0);
end
