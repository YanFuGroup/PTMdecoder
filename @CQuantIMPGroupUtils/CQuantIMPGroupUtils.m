classdef CQuantIMPGroupUtils
    % CQuantIMPGroupUtils
    % Utilities for IMP group quantification workflows.
    %
    % This class hosts mid-level quantification logic that composes
    % chromatogram utilities into higher-level IMP-group behaviors.
    methods (Static)
        
        esti_ratio = filter_and_normalize_peak_ratios(rt_grid, smoothed_intensity, esti_ratio, XIC_peaks, resFilterThres)

        [imp_max_props, peak_fwhms, ratio_each_XIC_peak, rt_bound] = ...
            compute_peak_features(rt_grid, smoothed_intensity, esti_ratio, XIC_peaks)

        idx_selected = select_best_peak_per_imp(imp_max_props, ratio_each_XIC_peak)
    end
end
