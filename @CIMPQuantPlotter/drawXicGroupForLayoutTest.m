function diagnostics = drawXicGroupForLayoutTest(ric, total_xic, categorized_intervals, ...
    current_imp_name, colors, file_base_path, layout)
% Draw a synthetic XIC group for layout tests.
color_map = containers.Map();
legend_map = containers.Map();
for idx = 1:numel(current_imp_name)
    color_map(current_imp_name{idx}) = colors(idx, :);
    legend_map(current_imp_name{idx}) = current_imp_name{idx};
end
diagnostics = plotXicGroupWithLayout(ric, total_xic, categorized_intervals, ...
    current_imp_name, file_base_path, color_map, legend_map, layout);
end
