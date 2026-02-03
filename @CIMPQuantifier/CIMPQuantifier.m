classdef CIMPQuantifier
    % Static utilities for IMP quantification

    methods (Static)
        % Quantify each group
        [has_nonzero_imp, imp_idx_nonzero, area_imp_final, rt_bound, idx_selected, ratio_each_XIC_peak] = ...
            quantGroup(cMs12DatasetIO, raw_name, ratio_raw, rt_raw, intensity_raw, low_mz_bound, high_mz_bound, selected_charge, minMSMSnum, alpha, resFilterThres)

        % Re-quantify each group
        [has_nonzero_imp, imp_idx_nonzero, area_imp_final, rt_bound, max_label, ratio_each_XIC_peak] = ...
            requantGroup(cMs12DatasetIO, raw_name, ratio_raw, rt_raw, intensity_raw, low_mz_bound, high_mz_bound, selected_charge, current_imp_rt_range, minMSMSnum)
    end
end