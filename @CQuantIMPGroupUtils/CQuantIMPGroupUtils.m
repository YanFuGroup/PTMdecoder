classdef CQuantIMPGroupUtils
    % CQuantIMPGroupUtils
    % Utilities for IMP group quantification workflows.
    %
    % This class hosts mid-level quantification logic that composes
    % chromatogram utilities into higher-level IMP-group behaviors.
    methods (Static)
        
        esti_ratio = filter_and_normalize_peak_ratios(rt_grid, smoothed_intensity, esti_ratio, XIC_peaks, resFilterThres)

        [imp_max_props, peak_fwhms, area_each_XIC_peak, rt_bound] = ...
            compute_peak_features(rt_grid, smoothed_intensity, esti_ratio, XIC_peaks)

        idx_selected = select_best_peak_per_imp(imp_max_props, area_each_XIC_peak)

        esti_ratio = refine_ratios_by_selection(esti_ratio, XIC_peaks, idx_selected)

        auxic = compute_final_area(rt_grid, smoothed_intensity, esti_ratio, XIC_peaks, idx_selected)

        auxic = compute_final_area_from_peak_areas(area_each_XIC_peak, idx_selected)

        [auxic, rt_bound, ratio_each_XIC_peak] = compute_imp_peak_area_and_ratio(...
            rt_grid, smoothed_intensity, esti_ratio, peak_ranges, final_XIC_peak_for_IMP, is_skip_vec)

        ric = build_ric_from_peaks(rt_grid, smoothed_intensity, esti_ratio, peak_ranges, is_skip_vec)

        [sort_rts, sort_ratioMatrix, rt_grid, smoothed_intensity, intensity, is_valid] = ...
            prepare_ms1_xic(cMs12DatasetIO, raw_name, current_rts, current_inten, current_ratioMatrix, ...
                minMSMSnum, low_mz_bound, high_mz_bound, selected_charge)

        [final_XIC_peak_for_IMP, max_label, is_skip_vec, peak_ranges] = ...
            prepare_peak_ranges_from_iso_rt_range(rt_grid, current_iso_rt_range, rt_error_tol)

        [idxNonZero, auxic, rt_bound, varargout] = filter_nonzero_xic(auxic, rt_bound, varargin)
    end
end
