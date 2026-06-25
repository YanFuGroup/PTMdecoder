function margin_px = measureRenderedLegendRightMarginPx(f, legend_column_start, figure_width_px)
% Measure right-side whitespace of rendered legend ink from the current figure frame.
frame = getframe(f);
gray_image = rgb2gray(frame.cdata);
img_h = size(gray_image, 1);
img_w = size(gray_image, 2);

if nargin < 3 || isempty(figure_width_px)
    figure_width_px = img_w;
end

scale_x = img_w / figure_width_px;
legend_column_start = max(1, min(img_w, round(legend_column_start * scale_x)));
row_start = max(1, round(img_h * 0.15));
row_end = min(img_h, round(img_h * 0.85));
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

margin_px = margin_px / scale_x;
end
