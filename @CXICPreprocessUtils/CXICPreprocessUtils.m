classdef CXICPreprocessUtils
    % CXICPreprocessUtils
    % Preprocessing utilities for IMP quant workflows.
    methods (Static)
        [rt_sorted, ratio_sorted, xic_rt, xic_intensity_smoothed, xic_intensity_raw, is_valid] = ...
            prepare_ms1_xic(cMs12DatasetIO, raw_name, rt_raw, intensity_raw, ratio_raw, ...
                minMSMSnum, low_mz_bound, high_mz_bound, selected_charge)

        [max_label, is_skip_vec, xic_peak_idx_bounds] = ...
            prepare_peak_ranges_from_imp_rt_range(xic_rt, current_imp_rt_range, rt_error_tol)

        [xic_peak_idx_bounds, xic_peak_rt_bounds, is_ok] = ...
            build_peak_bounds_from_candidates(xic_rt, candidate_rt_peaks)
    end
end
