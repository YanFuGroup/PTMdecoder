function area_imp_final = compute_final_area(xic_rt, xic_intensity_smoothed, ratio_estimated, XIC_peaks, idx_selected)
% compute_final_area
% Compute final area per IMP within the selected peak region.
% Deprecated: use compute_final_area_from_peak_areas instead.
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
%   idx_selected (K x 1 double)
%       Selected peak index per IMP
%
% Output:
%   area_imp_final (K x 1 double) area
%       Final area per IMP

num_imp = size(ratio_estimated, 2);
intensityMatrix = ratio_estimated .* xic_intensity_smoothed;
area_imp_final = zeros(num_imp, 1);

for idx_imp = 1:num_imp
    sel_peak_idx = idx_selected(idx_imp);
    idx_start = XIC_peaks(sel_peak_idx).left_bound;
    idx_end = XIC_peaks(sel_peak_idx).right_bound;
    area_imp_final(idx_imp, 1) = CChromatogramUtils.calculate_area(...
        xic_rt, intensityMatrix(:, idx_imp), idx_start, idx_end);
end
end
