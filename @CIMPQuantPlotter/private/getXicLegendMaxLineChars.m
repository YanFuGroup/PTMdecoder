function max_line_chars = getXicLegendMaxLineChars(layout)
% Estimate maximum legend label line length from layout pixel budget.
if isfield(layout, 'legend_default_icon_width_px') && ~isempty(layout.legend_default_icon_width_px)
    icon_width_px = layout.legend_default_icon_width_px;
else
    icon_width_px = 30;
end
text_width_px = layout.legend_max_width_px - layout.legend_padding_px - icon_width_px;
text_width_px = max(1, text_width_px);
if isfield(layout, 'legend_chars_per_px') && ~isempty(layout.legend_chars_per_px)
    chars_per_px = layout.legend_chars_per_px;
else
    chars_per_px = 0.15;
end
max_line_chars = floor(text_width_px * chars_per_px);
if isfield(layout, 'legend_min_line_chars') && ~isempty(layout.legend_min_line_chars)
    max_line_chars = max(layout.legend_min_line_chars, max_line_chars);
else
    max_line_chars = max(1, max_line_chars);
end
end
