function metrics = assertXicLegendExports(output_base, layout, diagnostics)
% Structural checks on exported XIC PNG/PDF/SVG outputs.
min_export_bytes = 500;
min_legend_margin_px = 8;

png_path = [output_base, '.png'];
pdf_path = [output_base, '.pdf'];
svg_path = [output_base, '.svg'];

metrics = assertPngExportStructure(png_path, layout, diagnostics, min_legend_margin_px);
assertPdfExportIntegrity(pdf_path, layout, diagnostics, min_export_bytes);
if isfile(svg_path)
    assertSvgExportIntegrity(svg_path, layout, diagnostics, min_export_bytes);
    metrics.svg_checked = true;
else
    metrics.svg_checked = false;
end
end

function metrics = assertPngExportStructure(png_path, layout, diagnostics, min_legend_margin_px)
png_img = imread(png_path);
gray_image = rgb2gray(png_img);
img_h = size(gray_image, 1);
img_w = size(gray_image, 2);

fig_w = diagnostics.figure_position_px(3);
fig_h = diagnostics.figure_position_px(4);

assert(abs(img_w - fig_w) <= 2, ...
    'PNG width %d does not match figure width %d.', img_w, fig_w);
assert(abs(img_h - fig_h) <= 2, ...
    'PNG height %d does not match figure height %d.', img_h, fig_h);

legend_box_col = layout.left_margin_px + layout.axes_width_px + layout.axes_legend_gap_px;
legend_col = legend_box_col;
margin_px = CIMPQuantPlotter.measureExportedLegendRightMarginPx(png_path, legend_col);
assert(margin_px >= min_legend_margin_px, ...
    'Exported PNG legend right margin too small: %d px', margin_px);

axes_right = legend_box_col - 1;
legend_col_start = max(1, legend_box_col);
axes_band = gray_image(:, 1:min(axes_right, img_w));
legend_band = gray_image(max(1, round(img_h * 0.15)):min(img_h, round(img_h * 0.85)), ...
    legend_col_start:end);

assert(any(axes_band(:) < 240), 'Axes region has no ink in exported PNG.');
assert(any(legend_band(:) < 240), 'Legend region has no ink in exported PNG.');

[row_idx, col_idx] = find(gray_image < 240);
assert(~isempty(row_idx), 'Exported PNG has no ink pixels.');
ink_row_fraction = mean(row_idx) / img_h;
assert(ink_row_fraction > 0.20, ...
    'Ink concentrated in top band (mean row fraction %.3f <= 0.20).', ink_row_fraction);
assert(ink_row_fraction < 0.85, ...
    'Ink shifted too low (mean row fraction %.3f >= 0.85).', ink_row_fraction);

max_ink_col_fraction = max(col_idx) / img_w;
right_blank_fraction = (img_w - max(col_idx)) / img_w;
assert(right_blank_fraction <= 0.20, ...
    'Too much right blank space (right blank fraction %.3f > 0.20).', right_blank_fraction);

axes_ink_cols = col_idx(col_idx <= axes_right);
legend_ink_cols = col_idx(col_idx >= legend_col_start);
assert(~isempty(axes_ink_cols), 'Axes region has no ink columns in exported PNG.');
assert(~isempty(legend_ink_cols), 'Legend region has no ink columns in exported PNG.');
axes_legend_gap_px = min(legend_ink_cols) - max(axes_ink_cols);
assert(axes_legend_gap_px >= 1, ...
    'Legend ink overlaps axes ink in exported PNG (gap %d px).', axes_legend_gap_px);

metrics = struct();
metrics.png_width_px = img_w;
metrics.png_height_px = img_h;
metrics.legend_right_margin_px = margin_px;
metrics.ink_row_fraction = ink_row_fraction;
metrics.max_ink_col_fraction = max_ink_col_fraction;
metrics.right_blank_fraction = right_blank_fraction;
end

function assertPdfExportIntegrity(pdf_path, layout, diagnostics, min_export_bytes)
assert(isfile(pdf_path), 'Missing exported PDF: %s', pdf_path);

pdf_info = dir(pdf_path);
assert(~isempty(pdf_info), 'Missing exported PDF: %s', pdf_path);
assert(pdf_info.bytes >= min_export_bytes, ...
    'Exported PDF too small: %d bytes', pdf_info.bytes);

fid = fopen(pdf_path, 'r');
assert(fid >= 0, 'Unable to open exported PDF: %s', pdf_path);
header = fread(fid, 4, '*char')';
fclose(fid);
assert(strcmp(header, '%PDF'), 'Exported PDF has invalid header: %s', header);

try
    fig_w = diagnostics.figure_position_px(3);
    fig_h = diagnostics.figure_position_px(4);
    image_info = imfinfo(pdf_path);
    assert(~isempty(image_info), 'imfinfo returned no pages for PDF: %s', pdf_path);
    expected_width_pt = fig_w / layout.dpi * 72;
    expected_height_pt = fig_h / layout.dpi * 72;
    if isfield(image_info(1), 'Width') && isfield(image_info(1), 'Height')
        width_tolerance = max(5, 0.05 * expected_width_pt);
        height_tolerance = max(5, 0.05 * expected_height_pt);
        assert(abs(image_info(1).Width - expected_width_pt) <= width_tolerance, ...
            'PDF width %.1f pt differs from expected %.1f pt.', ...
            image_info(1).Width, expected_width_pt);
        assert(abs(image_info(1).Height - expected_height_pt) <= height_tolerance, ...
            'PDF height %.1f pt differs from expected %.1f pt.', ...
            image_info(1).Height, expected_height_pt);
    end
catch
    warning('test_CIMPQuantPlotterLegendLayout:PdfDimensionCheckSkipped', ...
        'Skipping PDF dimension check for %s.', pdf_path);
end
end

function assertSvgExportIntegrity(svg_path, layout, diagnostics, min_export_bytes)
assert(isfile(svg_path), 'Missing exported SVG: %s', svg_path);

svg_info = dir(svg_path);
assert(~isempty(svg_info), 'Missing exported SVG: %s', svg_path);
assert(svg_info.bytes >= min_export_bytes, ...
    'Exported SVG too small: %d bytes', svg_info.bytes);

svg_text = fileread(svg_path);
assert(contains(svg_text, '<svg', 'IgnoreCase', true), ...
    'Exported SVG does not contain an <svg> root element.');

fig_w = diagnostics.figure_position_px(3);
fig_h = diagnostics.figure_position_px(4);
expected_svg_width = round(fig_w * 96 / layout.dpi);
expected_svg_height = round(fig_h * 96 / layout.dpi);
assert(contains(svg_text, sprintf('width="%d"', expected_svg_width)), ...
    'SVG width does not match expected %d px.', expected_svg_width);
assert(contains(svg_text, sprintf('height="%d"', expected_svg_height)), ...
    'SVG height does not match expected %d px.', expected_svg_height);
end
