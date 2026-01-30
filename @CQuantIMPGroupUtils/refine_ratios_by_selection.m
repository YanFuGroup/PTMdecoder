function esti_ratio = refine_ratios_by_selection(esti_ratio, XIC_peaks, idx_selected)
% refine_ratios_by_selection
% Keep only the selected peak region for each IMP and zero out others.
%
% Inputs:
%   esti_ratio   Estimated ratio matrix (len(rt_grid) x num_imp)
%   XIC_peaks    Struct array with left_bound/right_bound indices
%   idx_selected Selected peak index per IMP
%
% Output:
%   esti_ratio   Refined ratio matrix

num_imp = size(esti_ratio, 2);
keep_mask = false(size(esti_ratio));

for idx_imp = 1:num_imp
    sel_peak_idx = idx_selected(idx_imp);
    if sel_peak_idx > 0 && sel_peak_idx <= length(XIC_peaks)
        p_start = XIC_peaks(sel_peak_idx).left_bound;
        p_end   = XIC_peaks(sel_peak_idx).right_bound;
        keep_mask(p_start:p_end, idx_imp) = true;
    end
end

esti_ratio(~keep_mask) = 0;
end
