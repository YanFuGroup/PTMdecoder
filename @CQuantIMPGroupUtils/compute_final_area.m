function auxic = compute_final_area(rt_grid, smoothed_intensity, esti_ratio, XIC_peaks, idx_selected)
% compute_final_area
% Compute final area per IMP within the selected peak region.
% Deprecated: use compute_final_area_from_peak_areas instead.
%
% Inputs:
%   rt_grid (N x 1 double) minutes
%       RT grid vector
%   smoothed_intensity (N x 1 double) intensity
%       Smoothed XIC intensity (aligned to rt_grid)
%   esti_ratio (N x K double)
%       Estimated ratio matrix for K IMPs
%   XIC_peaks (1 x P struct)
%       Struct array with fields: left_bound/right_bound (indices into rt_grid)
%   idx_selected (K x 1 double)
%       Selected peak index per IMP
%
% Output:
%   auxic (K x 1 double) area
%       Final area per IMP

num_imp = size(esti_ratio, 2);
intensityMatrix = esti_ratio .* smoothed_intensity;
auxic = zeros(num_imp, 1);

for idx_imp = 1:num_imp
    sel_peak_idx = idx_selected(idx_imp);
    idx_start = XIC_peaks(sel_peak_idx).left_bound;
    idx_end = XIC_peaks(sel_peak_idx).right_bound;
    auxic(idx_imp, 1) = CChromatogramUtils.calculate_area(...
        rt_grid, intensityMatrix(:, idx_imp), idx_start, idx_end);
end
end
