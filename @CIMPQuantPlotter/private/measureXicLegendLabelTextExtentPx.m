function text_extent = measureXicLegendLabelTextExtentPx(ax, label, font_size, font_name)
% Measure a single legend label line extent in axes pixel coordinates.
text_handle = text(ax, 0, 0, label, 'Visible', 'off', 'FontSize', font_size, ...
    'FontName', font_name, 'Interpreter', 'none');
drawnow;
set(text_handle, 'Units', 'pixels');
text_extent = get(text_handle, 'Extent');
delete(text_handle);
end
