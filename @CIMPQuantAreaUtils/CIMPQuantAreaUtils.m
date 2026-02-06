classdef CIMPQuantAreaUtils
    % CIMPQuantAreaUtils
    % Area/aggregation utilities for IMP quant workflows.
    methods (Static)
        area_imp_final = compute_final_area(xic_rt, xic_intensity_smoothed, xic_ratio_estimated, xic_peak_idx_bounds, idx_selected)

        area_imp_final = get_final_area_from_peak_areas(area_imp_by_peak, idx_selected)

        [area_imp_final, xic_peak_rt_bounds, ratio_each_XIC_peak] = compute_imp_peak_area_and_ratio(...
            xic_rt, xic_intensity_smoothed, xic_ratio_estimated, xic_peak_idx_bounds, final_xic_peak_rt_bounds, is_skip_vec)

        ric = build_ric_from_peaks(xic_rt, xic_intensity_smoothed, xic_ratio_estimated, xic_peak_idx_bounds, is_skip_vec)

        [imp_idx_nonzero, area_imp_final, xic_peak_rt_bounds, varargout] = filter_nonzero_xic(area_imp_final, xic_peak_rt_bounds, varargin)
    end
end
