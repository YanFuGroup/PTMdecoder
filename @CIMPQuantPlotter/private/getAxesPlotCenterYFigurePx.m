function center_y = getAxesPlotCenterYFigurePx(ax, ax_pos)
% Return the vertical center of the axes plot box in figure pixel coordinates.
saved_units = get(ax, 'Units');
set(ax, 'Units', 'normalized');
inner_pos = get(ax, 'InnerPosition');
set(ax, 'Units', saved_units);

inner_bottom = ax_pos(2) + inner_pos(2) * ax_pos(4);
inner_height = inner_pos(4) * ax_pos(4);
center_y = inner_bottom + inner_height / 2;
end
