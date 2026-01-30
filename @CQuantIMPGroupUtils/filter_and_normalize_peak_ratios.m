function esti_ratio = filter_and_normalize_peak_ratios(rt_grid, smoothed_intensity, esti_ratio, XIC_peaks, resFilterThres)
% filter_and_normalize_peak_ratios
% For each detected XIC peak, remove IMPs with very small area
% and normalize ratios within the peak.
%
% Inputs:
%   rt_grid            RT grid vector
%   smoothed_intensity Smoothed XIC intensity (column vector with len(rt_grid) rows)
%   esti_ratio         Estimated ratio matrix (len(rt_grid) x num_imp)
%   XIC_peaks          Struct array with left_bound/right_bound indices
%   resFilterThres     Threshold (relative to max area in a peak)
%
% Output:
%   esti_ratio         Updated ratio matrix after filtering/normalization

num_imp = size(esti_ratio, 2);
intensityMatrix = esti_ratio .* smoothed_intensity;

for i_Xp = 1:length(XIC_peaks)
    curr_start = XIC_peaks(i_Xp).left_bound;
    curr_end = XIC_peaks(i_Xp).right_bound;

    % Calculate area for each IMP in this peak
    area_filter = zeros(num_imp, 1);
    for idx_imp = 1:num_imp
        area_filter(idx_imp) = CChromatogramUtils.calculate_area(...
            rt_grid, intensityMatrix(:, idx_imp), ...
            curr_start, curr_end);
    end
    
    % Filter: keep only IMPs with area >= max_area * threshold
    max_area = max(area_filter);
    keep_mask = area_filter >= max_area * resFilterThres;
    esti_ratio(curr_start:curr_end, ~keep_mask) = 0;

    % Normalize rows
    row_sum = sum(esti_ratio(curr_start:curr_end, :), 2);
    row_sum(row_sum == 0) = 1;
    esti_ratio(curr_start:curr_end, :) = esti_ratio(curr_start:curr_end, :) ./ repmat(row_sum, 1, num_imp);
end
end
