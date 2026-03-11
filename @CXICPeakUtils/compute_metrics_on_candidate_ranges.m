function [imp_max_props, area_imp_by_peak, ratio_each_XIC_peak] = compute_metrics_on_candidate_ranges(...
    xic_rt, xic_intensity_smoothed, rt_sorted, ratio_sorted, xic_peak_idx_bounds, resFilterThres)
% Compute candidate-peak metrics from kernel-estimated ratios.
% Input:
%   xic_rt/xic_intensity_smoothed
%       XIC grid and smoothed intensity
%   rt_sorted/ratio_sorted
%       PSM RTs and ratio matrix used for kernel estimation
%   xic_peak_idx_bounds
%       Candidate peak index intervals
%   resFilterThres
%       Relative area filter threshold used during ratio filtering
% Output:
%   imp_max_props (K x P double)
%       Max ratio_estimated per IMP per candidate interval
%   area_imp_by_peak (K x P double)
%       IMP-specific area under (ratio_estimated .* XIC) per interval
%   ratio_each_XIC_peak (K x P double)
%       IMP area divided by total XIC area for each candidate interval

num_imp = size(ratio_sorted, 2);
num_peaks = numel(xic_peak_idx_bounds);

imp_max_props = zeros(num_imp, num_peaks);
area_imp_by_peak = zeros(num_imp, num_peaks);
ratio_each_XIC_peak = zeros(num_imp, num_peaks);

xic_ratio_estimated = CXICPeakUtils.calculate_kernel_ratio(...
    xic_rt, rt_sorted, ratio_sorted, xic_peak_idx_bounds, true);
xic_ratio_estimated = CXICPeakUtils.filter_and_normalize_peak_ratios(...
    xic_rt, xic_intensity_smoothed, xic_ratio_estimated, xic_peak_idx_bounds, resFilterThres);

intensity_matrix = xic_ratio_estimated .* xic_intensity_smoothed;

for idx_peak = 1:num_peaks
    idx_start = xic_peak_idx_bounds(idx_peak).idx_start;
    idx_end = xic_peak_idx_bounds(idx_peak).idx_end;
    peak_total_area = CXICSignalUtils.calculate_area(...
        xic_rt, xic_intensity_smoothed, idx_start, idx_end);
    if peak_total_area <= 0
        continue;
    end

    ratio_slice = xic_ratio_estimated(idx_start:idx_end, :);
    imp_max_props(:, idx_peak) = max(ratio_slice, [], 1)';

    for idx_imp = 1:num_imp
        area_imp_by_peak(idx_imp, idx_peak) = CXICSignalUtils.calculate_area(...
            xic_rt, intensity_matrix(:, idx_imp), idx_start, idx_end);
        ratio_each_XIC_peak(idx_imp, idx_peak) = ...
            area_imp_by_peak(idx_imp, idx_peak) / peak_total_area;
    end
end
end