classdef CChromatogramUtils
    % CChromatogramUtils
    % A utility class for chromatogram signal processing, including
    % smoothing, peak detection, and data preprocessing.
    
    methods (Static)
        [sort_rts, sort_ratioMatrix, is_valid] = preprocess_ms1_inputs(current_rts, current_inten, current_ratioMatrix, minMSMSnum)
        
        [rt_grid, smoothed_intensity, intensity] = get_smoothed_xic(ms12DatasetIO, raw_name, low_mz_bound, high_mz_bound, selected_charge)
        
        fwhm = get_fwhm(peak_rts, peak_intens)
        
        XIC_peaks = detect_xic_peaks(rt_grid, smoothed_intensity, raw_intensity, sort_rts, alpha)
        
        esti_ratio = calculate_kernel_ratio(rt_grid, sort_rts, sort_ratioMatrix, peak_ranges, is_broadcast)
        
        [final_XIC_peak, max_label, is_skip_vec] = parse_imp_rt_ranges(imp_rt_range, is_skip_vec)
        
        peak_ranges = map_rt_to_indices(rt_grid, final_XIC_peak, is_skip_vec, rt_error_tol)

        [rec_rt, rec_inten] = get_closed_peak_data(rt_grid, intensity_full, idx_start, idx_end)
        
        area = calculate_area(rt_grid, intensity_full, idx_start, idx_end)
    end
end
