classdef CQuantIMPGroupPeakUtils
    % CQuantIMPGroupPeakUtils
    % Peak detection/selection utilities for IMP group quant workflows.
    methods (Static)
        esti_ratio = filter_and_normalize_peak_ratios(rt_grid, smoothed_intensity, esti_ratio, XIC_peaks, resFilterThres)

        [imp_max_props, peak_fwhms, area_each_XIC_peak, rt_bound] = ...
            compute_peak_features(rt_grid, smoothed_intensity, esti_ratio, XIC_peaks)

        idx_selected = select_best_peak_per_imp(imp_max_props, area_each_XIC_peak)

        esti_ratio = refine_ratios_by_selection(esti_ratio, XIC_peaks, idx_selected)
    end
end
