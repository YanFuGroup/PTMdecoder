classdef CQuantIMPGroupAreaUtils
    % CQuantIMPGroupAreaUtils
    % Area/aggregation utilities for IMP group quant workflows.
    methods (Static)
        area_imp_final = compute_final_area(xic_rt, xic_intensity_smoothed, ratio_estimated, XIC_peaks, idx_selected)

        area_imp_final = get_final_area_from_peak_areas(area_imp_by_peak, idx_selected)

        [area_imp_final, rt_bound, ratio_each_XIC_peak] = compute_imp_peak_area_and_ratio(...
            xic_rt, xic_intensity_smoothed, ratio_estimated, peak_ranges, final_XIC_peak_for_IMP, is_skip_vec)

        ric = build_ric_from_peaks(xic_rt, xic_intensity_smoothed, ratio_estimated, peak_ranges, is_skip_vec)

        [imp_idx_nonzero, area_imp_final, rt_bound, varargout] = filter_nonzero_xic(area_imp_final, rt_bound, varargin)
    end
end
