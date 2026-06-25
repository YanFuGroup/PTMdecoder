function exportXicFigureWithLayout(f, file_base_path, layout, fig_w, fig_h)
% Export XIC figure to PNG, PDF, and SVG at logical layout dimensions.
% PNG pixels match figure pixels (PaperSize in inches = fig/dpi, print -r dpi).
if nargin < 5 || isempty(fig_w) || isempty(fig_h)
    fig_pos = get(f, 'Position');
    fig_w = fig_pos(3);
    fig_h = fig_pos(4);
end

png_file = [file_base_path, '.png'];
pdf_file = [file_base_path, '.pdf'];
svg_file = [file_base_path, '.svg'];

set(f, 'Units', 'pixels', 'Position', [50, 50, fig_w, fig_h]);
set(f, 'PaperUnits', 'inches', ...
    'PaperPosition', [0, 0, fig_w / layout.dpi, fig_h / layout.dpi], ...
    'PaperSize', [fig_w / layout.dpi, fig_h / layout.dpi], ...
    'PaperPositionMode', 'manual', 'InvertHardcopy', 'off', 'Renderer', 'opengl');
drawnow;
print(f, png_file, '-dpng', sprintf('-r%d', layout.dpi));

set(f, 'Renderer', 'painters');
drawnow;
print(f, pdf_file, '-dpdf', '-painters');
print(f, svg_file, '-dsvg', '-painters');
end
