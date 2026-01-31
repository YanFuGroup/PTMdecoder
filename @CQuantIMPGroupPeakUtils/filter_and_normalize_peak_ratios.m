function ratio_estimated = filter_and_normalize_peak_ratios(xic_rt, xic_intensity_smoothed, ratio_estimated, XIC_peaks, resFilterThres)
% filter_and_normalize_peak_ratios
% For each detected XIC peak, remove IMPs with very small area
% and normalize ratios within the peak.
%
% Inputs:
%   xic_rt (N x 1 double) minutes
%       RT grid vector
%   xic_intensity_smoothed (N x 1 double) intensity
%       Smoothed XIC intensity (aligned to xic_rt)
%   ratio_estimated (N x K double)
%       Estimated ratio matrix for K IMPs
%   XIC_peaks (1 x P struct)
%       Struct array with fields: left_bound/right_bound (indices into xic_rt)
%   resFilterThres (1 x 1 double)
%       Threshold (relative to max area in a peak)
%
% Output:
%   ratio_estimated (N x K double)
%       Updated ratio matrix after filtering/normalization

num_imp = size(ratio_estimated, 2);
intensityMatrix = ratio_estimated .* xic_intensity_smoothed;

for i_Xp = 1:length(XIC_peaks)
    curr_start = XIC_peaks(i_Xp).left_bound;
    curr_end = XIC_peaks(i_Xp).right_bound;

    % Calculate area for each IMP in this peak
    area_filter = zeros(num_imp, 1);
    for idx_imp = 1:num_imp
        area_filter(idx_imp) = CChromatogramUtils.calculate_area(...
                xic_rt, intensityMatrix(:, idx_imp), ...
            curr_start, curr_end);
    end
    
    % Filter: keep only IMPs with area >= max_area * threshold
    max_area = max(area_filter);
    keep_mask = area_filter >= max_area * resFilterThres;
    ratio_estimated(curr_start:curr_end, ~keep_mask) = 0;

    % Normalize rows
    row_sum = sum(ratio_estimated(curr_start:curr_end, :), 2);
    row_sum(row_sum == 0) = 1;
    ratio_estimated(curr_start:curr_end, :) = ratio_estimated(curr_start:curr_end, :) ./ repmat(row_sum, 1, num_imp);
end
end
