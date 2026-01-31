classdef CQuantIMPGroupPreprocessUtils
    % CQuantIMPGroupPreprocessUtils
    % Preprocessing utilities for IMP group quant workflows.
    methods (Static)
        [sort_rts, sort_ratioMatrix, rt_grid, smoothed_intensity, intensity, is_valid] = ...
            prepare_ms1_xic(cMs12DatasetIO, raw_name, current_rts, current_inten, current_ratioMatrix, ...
                minMSMSnum, low_mz_bound, high_mz_bound, selected_charge)

        [final_XIC_peak_for_IMP, max_label, is_skip_vec, peak_ranges] = ...
            prepare_peak_ranges_from_imp_rt_range(rt_grid, current_iso_rt_range, rt_error_tol)
    end
end
