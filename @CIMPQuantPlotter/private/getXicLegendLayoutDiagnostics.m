function diagnostics = getXicLegendLayoutDiagnostics(f, ax, h_legend, layout)
% Inspect rendered legend geometry in figure pixel coordinates.
set(f, 'Units', 'pixels');
set(ax, 'Units', 'pixels');
drawnow;

fig_pos = get(f, 'Position');
ax_pos = get(ax, 'Position');
fig_w = fig_pos(3);
fig_h = fig_pos(4);
axes_right_with_gap = ax_pos(1) + ax_pos(3) + layout.axes_legend_gap_px;

text_bounds = measureXicLegendRenderedBoundsFigurePx(ax, h_legend);
ink_bounds = measureXicLegendInkBoundsFigurePx(f, ax, layout, fig_w, fig_h);
print_width_scale = 1.0;
if isfield(layout, 'legend_print_width_scale') && ~isempty(layout.legend_print_width_scale)
    print_width_scale = layout.legend_print_width_scale;
end
legend_pos = text_bounds.legend_position_px;
legend_right = legend_pos(1) + legend_pos(3);
legend_top = legend_pos(2) + legend_pos(4);
text_extents = text_bounds.text_extents;

min_text_left = text_bounds.min_left;
max_text_right = text_bounds.max_right;
min_text_bottom = text_bounds.min_bottom;
max_text_top = text_bounds.max_top;
if ~isempty(text_extents) && any(text_extents(:, 3) > 0)
    text_left_edge = min(text_extents(:, 1));
else
    text_left_edge = legend_pos(1) + text_bounds.icon_width;
end

diagnostics = struct();
diagnostics.figure_position_px = fig_pos;
diagnostics.axes_position_px = ax_pos;
diagnostics.legend_position_px = legend_pos;
diagnostics.legend_text_extents_px = text_extents;
diagnostics.legend_center_y_px = ink_bounds.center_y;
diagnostics.axes_center_y_px = getAxesPlotCenterYFigurePx(ax, ax_pos);
diagnostics.is_legend_inside_figure = legend_pos(1) >= 0 && legend_pos(2) >= 0 && ...
    legend_right <= fig_w && legend_top <= fig_h;
diagnostics.is_text_inside_figure = min_text_left >= 0 && min_text_bottom >= 0 && ...
    max_text_right <= fig_w && max_text_top <= fig_h;
diagnostics.is_text_clear_of_axes = text_left_edge >= axes_right_with_gap;
diagnostics.max_text_right_px = max_text_right;
scaled_max_right = legend_pos(1) + (max_text_right - legend_pos(1)) * print_width_scale;
diagnostics.min_text_right_margin_px = fig_w - scaled_max_right;
diagnostics.summary = sprintf(['legend_fig=[%.1f %.1f %.1f %.1f], fig=[%.1f %.1f], ', ...
    'text_right=%.1f, text_left=%.1f, axes_right_gap=%.1f, margin=%.1f, center_y=%.1f'], ...
    legend_pos(1), legend_pos(2), legend_pos(3), legend_pos(4), ...
    fig_w, fig_h, max_text_right, text_left_edge, axes_right_with_gap, ...
    diagnostics.min_text_right_margin_px, ink_bounds.center_y);
end
