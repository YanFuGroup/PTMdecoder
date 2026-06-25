function assertXicExportIntegrity(file_base_path, layout, fig_w, fig_h)
% Validate exported XIC figure files after print.
% Throws CIMPQuantPlotter:LegendLayoutClipped for PNG layout issues and
% CIMPQuantPlotter:LegendExportInvalid for corrupt PDF/SVG output.
min_export_bytes = 500;
min_legend_margin_px = 8;

png_file = [file_base_path, '.png'];
pdf_file = [file_base_path, '.pdf'];
svg_file = [file_base_path, '.svg'];

try
    assertPngExportIntegrity(png_file, layout, fig_w, fig_h, min_legend_margin_px);
catch export_error
    if strcmp(export_error.identifier, 'CIMPQuantPlotter:LegendExportInvalid')
        error('CIMPQuantPlotter:LegendLayoutClipped', '%s', export_error.message);
    end
    rethrow(export_error);
end

assertPdfExportIntegrity(pdf_file, layout, fig_w, fig_h, min_export_bytes);
assertSvgExportIntegrity(svg_file, layout, fig_w, fig_h, min_export_bytes);
end

function assertPngExportIntegrity(png_file, layout, fig_w, fig_h, min_legend_margin_px)
if ~isfile(png_file)
    error('CIMPQuantPlotter:LegendExportInvalid', 'Missing exported PNG: %s', png_file);
end

png_img = imread(png_file);
gray_image = rgb2gray(png_img);
img_h = size(gray_image, 1);
img_w = size(gray_image, 2);

if abs(img_w - fig_w) > 2
    error('CIMPQuantPlotter:LegendExportInvalid', ...
        'PNG width %d does not match figure width %d.', img_w, fig_w);
end
if abs(img_h - fig_h) > 2
    error('CIMPQuantPlotter:LegendExportInvalid', ...
        'PNG height %d does not match figure height %d.', img_h, fig_h);
end

legend_col = layout.left_margin_px + layout.axes_width_px + layout.axes_legend_gap_px;
margin_px = CIMPQuantPlotter.measureExportedLegendRightMarginPx(png_file, legend_col);
if margin_px < min_legend_margin_px
    error('CIMPQuantPlotter:LegendExportInvalid', ...
        'Exported PNG legend right margin too small: %d px', margin_px);
end

axes_right = legend_col - 1;
legend_col_start = max(1, legend_col);
axes_band = gray_image(:, 1:min(axes_right, img_w));
legend_band = gray_image(max(1, round(img_h * 0.15)):min(img_h, round(img_h * 0.85)), ...
    legend_col_start:end);
if ~any(axes_band(:) < 240)
    error('CIMPQuantPlotter:LegendExportInvalid', 'Axes region has no ink in exported PNG.');
end
if ~any(legend_band(:) < 240)
    error('CIMPQuantPlotter:LegendExportInvalid', 'Legend region has no ink in exported PNG.');
end
end

function assertPdfExportIntegrity(pdf_file, layout, fig_w, fig_h, min_export_bytes)
if ~isfile(pdf_file)
    error('CIMPQuantPlotter:LegendExportInvalid', 'Missing exported PDF: %s', pdf_file);
end

pdf_info = dir(pdf_file);
if pdf_info.bytes < min_export_bytes
    error('CIMPQuantPlotter:LegendExportInvalid', ...
        'Exported PDF too small: %d bytes', pdf_info.bytes);
end

fid = fopen(pdf_file, 'r');
if fid < 0
    error('CIMPQuantPlotter:LegendExportInvalid', 'Unable to open exported PDF: %s', pdf_file);
end
header = fread(fid, 4, '*char')';
fclose(fid);
if ~strcmp(header, '%PDF')
    error('CIMPQuantPlotter:LegendExportInvalid', ...
        'Exported PDF has invalid header: %s', header);
end

try
  image_info = imfinfo(pdf_file);
  if isempty(image_info)
      warning('CIMPQuantPlotter:LegendExportIntegrity', ...
          'imfinfo returned no pages for PDF: %s', pdf_file);
      return;
  end
  expected_width_pt = fig_w / layout.dpi * 72;
  expected_height_pt = fig_h / layout.dpi * 72;
  if isfield(image_info(1), 'Width') && isfield(image_info(1), 'Height')
      width_tolerance = max(5, 0.05 * expected_width_pt);
      height_tolerance = max(5, 0.05 * expected_height_pt);
      if abs(image_info(1).Width - expected_width_pt) > width_tolerance || ...
              abs(image_info(1).Height - expected_height_pt) > height_tolerance
          warning('CIMPQuantPlotter:LegendExportIntegrity', ...
              ['PDF page size [%.1f %.1f] differs from expected [%.1f %.1f] pt.'], ...
              image_info(1).Width, image_info(1).Height, expected_width_pt, expected_height_pt);
      end
  end
catch
  warning('CIMPQuantPlotter:LegendExportIntegrity', ...
      'Skipping PDF dimension check for %s.', pdf_file);
end
end

function assertSvgExportIntegrity(svg_file, layout, fig_w, fig_h, min_export_bytes)
if ~isfile(svg_file)
    error('CIMPQuantPlotter:LegendExportInvalid', 'Missing exported SVG: %s', svg_file);
end

svg_info = dir(svg_file);
if svg_info.bytes < min_export_bytes
    error('CIMPQuantPlotter:LegendExportInvalid', ...
        'Exported SVG too small: %d bytes', svg_info.bytes);
end

svg_text = fileread(svg_file);
if ~contains(svg_text, '<svg', 'IgnoreCase', true)
    error('CIMPQuantPlotter:LegendExportInvalid', ...
        'Exported SVG does not contain an <svg> root element.');
end

expected_svg_width = round(fig_w * 96 / layout.dpi);
expected_svg_height = round(fig_h * 96 / layout.dpi);
width_match = contains(svg_text, sprintf('width="%d"', expected_svg_width));
height_match = contains(svg_text, sprintf('height="%d"', expected_svg_height));
if ~(width_match && height_match)
    warning('CIMPQuantPlotter:LegendExportIntegrity', ...
        'SVG width/height [%d %d] differ from expected [%d %d] px.', ...
        extractSvgDimension(svg_text, 'width'), extractSvgDimension(svg_text, 'height'), ...
        expected_svg_width, expected_svg_height);
end
end

function dimension = extractSvgDimension(svg_text, attribute_name)
dimension = NaN;
pattern = sprintf('%s="(\\d+)"', attribute_name);
tokens = regexp(svg_text, pattern, 'tokens', 'once');
if ~isempty(tokens)
    dimension = str2double(tokens{1});
end
end
