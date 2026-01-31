classdef CQuantIMPGroupPreprocessUtils
    % CQuantIMPGroupPreprocessUtils
    % Preprocessing utilities for IMP group quant workflows.
    methods (Static)
        [rt_sorted, ratio_sorted, xic_rt, xic_intensity_smoothed, xic_intensity_raw, is_valid] = ...
            prepare_ms1_xic(cMs12DatasetIO, raw_name, rt_raw, intensity_raw, ratio_raw, ...
                minMSMSnum, low_mz_bound, high_mz_bound, selected_charge)

        [final_XIC_peak_for_IMP, max_label, is_skip_vec, peak_ranges] = ...
            prepare_peak_ranges_from_imp_rt_range(xic_rt, current_imp_rt_range, rt_error_tol)
    end
end
