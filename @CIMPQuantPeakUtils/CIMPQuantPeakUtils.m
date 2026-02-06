classdef CIMPQuantPeakUtils
    % CIMPQuantPeakUtils
    % Peak detection/selection utilities for IMP quant workflows.
    methods (Static)
        is_reserved = hasMinRows(ratio_matrix, min_rows)
        
        xic_ratio_estimated = filter_and_normalize_peak_ratios(xic_rt, xic_intensity_smoothed, xic_ratio_estimated, xic_peak_idx_bounds, resFilterThres)

        [imp_max_props, peak_fwhms, area_imp_by_peak, xic_peak_rt_bounds] = ...
            compute_peak_features(xic_rt, xic_intensity_smoothed, xic_ratio_estimated, xic_peak_idx_bounds)

            xic_ratio_estimated = calculate_kernel_ratio(xic_rt, rt_sorted, ratio_sorted, xic_peak_idx_bounds, is_broadcast)

        idx_selected = select_best_peak_per_imp(imp_max_props, area_imp_by_peak)

        xic_ratio_estimated = refine_ratios_by_selection(xic_ratio_estimated, xic_peak_idx_bounds, idx_selected)
    end
end
