classdef CQuantIMPGroupPeakUtils
    % CQuantIMPGroupPeakUtils
    % Peak detection/selection utilities for IMP group quant workflows.
    methods (Static)
        ratio_estimated = filter_and_normalize_peak_ratios(xic_rt, xic_intensity_smoothed, ratio_estimated, XIC_peaks, resFilterThres)

        [imp_max_props, peak_fwhms, area_imp_by_peak, rt_bound] = ...
            compute_peak_features(xic_rt, xic_intensity_smoothed, ratio_estimated, XIC_peaks)

        idx_selected = select_best_peak_per_imp(imp_max_props, area_imp_by_peak)

        ratio_estimated = refine_ratios_by_selection(ratio_estimated, XIC_peaks, idx_selected)
    end
end
