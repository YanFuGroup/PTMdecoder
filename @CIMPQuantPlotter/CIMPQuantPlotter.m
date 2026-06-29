classdef CIMPQuantPlotter
    % Plot XICs for IMP groups
    %
    % Static method prototypes below; implementations live in same-named .m files.

    methods (Static)
        drawGroup(ms12DatasetIO, minMSMSnum, raw_name, ratio_raw, rt_raw, ...
            intensity_raw, low_mz_bound, high_mz_bound, selected_charge, ...
            current_imp_rt_range, current_imp_name, dir_save, color_map, legend_map, xic_layout)
        layout = getXicLegendLayoutConfig()
        labels = prepareXicLegendLabels(labels, layout, max_line_chars)
        diagnostics = plotXicGroupWithLayout(ric, total_xic, categorized_intervals, ...
            current_imp_name, file_base_path, color_map, legend_map, layout)
    end

    methods (Hidden, Static)
        margin_px = measureExportedLegendRightMarginPx(png_path, legend_column_start)
    end
end
