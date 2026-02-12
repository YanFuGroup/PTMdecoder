classdef CXICAlignPeakSelector
    % Peak selection using alignment-aware scoring

    methods (Static)
        idx_selected = select_by_alignment(imp_max_props, area_imp_by_peak, xic_peak_rt_bounds, rt_pred, options)
    end
end
