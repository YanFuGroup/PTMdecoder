function xic_ratio_estimated = refine_ratios_by_selection(xic_ratio_estimated, xic_peak_idx_bounds, idx_selected)
% refine_ratios_by_selection
% Keep only the selected peak region for each IMP and zero out others.
% Deprecated: this function is no longer used in the main workflow since the final area is gotten from precomputed peak areas.
%
% Inputs:
%   xic_ratio_estimated (N x K double)
%       Estimated ratio matrix
%   xic_peak_idx_bounds (1 x P struct)
%       Struct array with idx_start/idx_end indices
%   idx_selected (K x 1 double)
%       Selected peak index per IMP
%
% Output:
%   xic_ratio_estimated (N x K double)
%       Refined ratio matrix

num_imp = size(xic_ratio_estimated, 2);
keep_mask = false(size(xic_ratio_estimated));

for idx_imp = 1:num_imp
    sel_peak_idx = idx_selected(idx_imp);
    if sel_peak_idx > 0 && sel_peak_idx <= length(xic_peak_idx_bounds)
        p_start = xic_peak_idx_bounds(sel_peak_idx).idx_start;
        p_end   = xic_peak_idx_bounds(sel_peak_idx).idx_end;
        keep_mask(p_start:p_end, idx_imp) = true;
    end
end

xic_ratio_estimated(~keep_mask) = 0;
end
