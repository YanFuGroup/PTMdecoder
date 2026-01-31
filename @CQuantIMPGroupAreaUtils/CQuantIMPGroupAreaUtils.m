classdef CQuantIMPGroupAreaUtils
    % CQuantIMPGroupAreaUtils
    % Area/aggregation utilities for IMP group quant workflows.
    methods (Static)
        auxic = compute_final_area(rt_grid, smoothed_intensity, esti_ratio, XIC_peaks, idx_selected)

        auxic = get_final_area_from_peak_areas(area_each_XIC_peak, idx_selected)

        [auxic, rt_bound, ratio_each_XIC_peak] = compute_imp_peak_area_and_ratio(...
            rt_grid, smoothed_intensity, esti_ratio, peak_ranges, final_XIC_peak_for_IMP, is_skip_vec)

        ric = build_ric_from_peaks(rt_grid, smoothed_intensity, esti_ratio, peak_ranges, is_skip_vec)

        [idxNonZero, auxic, rt_bound, varargout] = filter_nonzero_xic(auxic, rt_bound, varargin)
    end
end
