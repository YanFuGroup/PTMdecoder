function [color_map, legend_map] = buildXicLegendMaps(names, colors)
% Build color/legend maps for synthetic XIC layout tests.
color_map = containers.Map();
legend_map = containers.Map();
for idx = 1:numel(names)
    color_map(names{idx}) = colors(idx, :);
    legend_map(names{idx}) = names{idx};
end
end
