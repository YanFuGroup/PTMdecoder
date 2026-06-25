function bounds = measureXicLegendInkBoundsFigurePx(f, ax, layout, fig_w, fig_h)
% Measure legend vertical center from rendered figure ink pixels.
set(f, 'Units', 'pixels');
ax_pos = get(ax, 'Position');
if nargin < 5 || isempty(fig_h)
    fig_pos = get(f, 'Position');
    fig_w = fig_pos(3);
    fig_h = fig_pos(4);
end
legend_col_start = max(1, round(ax_pos(1) + ax_pos(3) + layout.axes_legend_gap_px * 0.25));

frame = getframe(f);
img = rgb2gray(frame.cdata);
img_h = size(img, 1);
img_w = size(img, 2);
size_tolerance_px = max(2, round(0.05 * max(fig_h, fig_w)));
if abs(img_h - fig_h) > size_tolerance_px || abs(img_w - fig_w) > size_tolerance_px
    bounds = struct('center_y', getAxesPlotCenterYFigurePx(ax, ax_pos));
    return;
end

scale_y = fig_h / img_h;
legend_band = img(:, legend_col_start:end);
[row_idx, ~] = find(legend_band < 240);
if isempty(row_idx)
    bounds = struct('center_y', getAxesPlotCenterYFigurePx(ax, ax_pos));
    return;
end

row_top = min(row_idx);
row_bottom = max(row_idx);
max_top = fig_h - (row_top - 1) * scale_y;
min_bottom = fig_h - (row_bottom - 1) * scale_y;
bounds = struct('center_y', (min_bottom + max_top) / 2);
end
