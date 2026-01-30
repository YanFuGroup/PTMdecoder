function auxic = compute_final_area(rt_grid, smoothed_intensity, esti_ratio, XIC_peaks, idx_selected)
% compute_final_area
% Compute final area per IMP within the selected peak region.
%
% Inputs:
%   rt_grid            RT grid vector
%   smoothed_intensity Smoothed XIC intensity (column vector)
%   esti_ratio         Estimated ratio matrix (len(rt_grid) x num_imp)
%   XIC_peaks          Struct array with left_bound/right_bound indices
%   idx_selected       Selected peak index per IMP
%
% Output:
%   auxic              Final area per IMP

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
