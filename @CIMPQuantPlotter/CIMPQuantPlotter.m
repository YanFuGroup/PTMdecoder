classdef CIMPQuantPlotter
    % Plot XICs for IMP groups

    methods (Static)
        drawGroup(ms12DatasetIO, minMSMSnum, raw_name, ratio_raw, rt_raw, ...
            intensity_raw, low_mz_bound, high_mz_bound, selected_charge, ...
            current_imp_rt_range, current_imp_name, dir_save, color_map, legend_map)
        layout = getXicLegendLayoutConfig()
        labels = prepareXicLegendLabels(labels, layout, max_line_chars)
        diagnostics = drawXicGroupForLayoutTest(ric, total_xic, categorized_intervals, ...
            current_imp_name, colors, file_base_path, layout)
    end
end
