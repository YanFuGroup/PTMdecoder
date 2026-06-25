function margin_px = measureExportedLegendRightMarginPx(png_path, legend_column_start)
% Measure the minimum right margin of legend ink in exported PNG pixels.
% legend_column_start is in exported image coordinates (logical layout * dpi/96).
image_data = imread(png_path);
gray_image = rgb2gray(image_data);
if nargin < 2 || isempty(legend_column_start)
    legend_column_start = 1;
end

img_h = size(gray_image, 1);
row_start = max(1, round(img_h * 0.15));
row_end = min(img_h, round(img_h * 0.85));
legend_column_start = max(1, min(legend_column_start, size(gray_image, 2)));
legend_band = gray_image(row_start:row_end, legend_column_start:end);
margin_px = size(legend_band, 2);
for row_idx = 1:size(legend_band, 1)
    row = legend_band(row_idx, :);
    ink_columns = find(row < 240);
    if isempty(ink_columns)
        continue;
    end
    row_margin = size(legend_band, 2) - max(ink_columns);
    margin_px = min(margin_px, row_margin);
end
end
