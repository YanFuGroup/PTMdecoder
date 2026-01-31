classdef CChromatogramUtils
    % CChromatogramUtils
    % A utility class for chromatogram signal processing, including
    % smoothing, peak detection, and data preprocessing.
    
    methods (Static)
        [rt_sorted, ratio_sorted, is_valid] = preprocess_ms1_inputs(rt_raw, intensity_raw, ratio_raw, minMSMSnum)
        
        [xic_rt, xic_intensity_smoothed, xic_intensity_raw] = get_smoothed_xic(ms12DatasetIO, raw_name, low_mz_bound, high_mz_bound, selected_charge)
        
        fwhm = get_fwhm(peak_rts, peak_intens)
        
        XIC_peaks = detect_xic_peaks(xic_rt, xic_intensity_smoothed, xic_intensity_raw, rt_sorted, alpha)
        
        ratio_estimated = calculate_kernel_ratio(xic_rt, rt_sorted, ratio_sorted, peak_ranges, is_broadcast)
        
        [final_XIC_peak, max_label, is_skip_vec] = parse_imp_rt_ranges(imp_rt_range, is_skip_vec)
        
        peak_ranges = map_rt_to_indices(xic_rt, final_XIC_peak, is_skip_vec, rt_error_tol)

        [rec_rt, rec_inten] = get_closed_peak_data(xic_rt, xic_intensity_full, idx_start, idx_end)
        
        area = calculate_area(xic_rt, xic_intensity_full, idx_start, idx_end)
    end
end
