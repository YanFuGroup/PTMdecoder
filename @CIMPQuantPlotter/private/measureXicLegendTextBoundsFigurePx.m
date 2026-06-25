function bounds = measureXicLegendTextBoundsFigurePx(ax, h_legend)
% Estimate rendered legend text bounds in figure pixel coordinates.
% Horizontal sizing uses proxy text metrics; vertical extents delegate to
% measureXicLegendRenderedBoundsFigurePx for top-down alignment.
set(ax, 'Units', 'pixels');
set(h_legend, 'Units', 'pixels');
legend_pos_ax = get(h_legend, 'Position');

rendered_bounds = measureXicLegendRenderedBoundsFigurePx(ax, h_legend);
strings = h_legend.String;
if ischar(strings)
    strings = cellstr(strings);
elseif isstring(strings)
    strings = cellstr(strings);
end

icon_width = rendered_bounds.icon_width;
font_size = h_legend.FontSize;
font_name = h_legend.FontName;
max_text_width = rendered_bounds.max_text_width;
for idx_string = 1:numel(strings)
    label_lines = splitXicLegendLabelLines(strings{idx_string});
    for idx_line = 1:numel(label_lines)
        text_extent = measureXicLegendLabelTextExtentPx( ...
            ax, label_lines{idx_line}, font_size, font_name);
        max_text_width = max(max_text_width, text_extent(3));
    end
end

auto_text_width = max(0, legend_pos_ax(3) - icon_width);
max_text_width = max(max_text_width, auto_text_width);

bounds = struct('icon_width', icon_width, 'max_text_width', max_text_width, ...
    'content_height', rendered_bounds.content_height, ...
    'max_right', rendered_bounds.max_right);
end
