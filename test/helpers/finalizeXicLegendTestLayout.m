function layout = finalizeXicLegendTestLayout(layout)
% Derive coherent figure_width_px, axes_width_px, and legend_max_width_px
% for shrunk test layouts.
%
% Tests can either provide axes_width_fraction with a fixed figure_width_px,
% or use the legacy explicit axes/legend widths when exercising edge cases.

if ~isfield(layout, 'legend_print_safety_px') || isempty(layout.legend_print_safety_px)
    layout.legend_print_safety_px = 32;
end

if isfield(layout, 'figure_width_px') && ~isempty(layout.figure_width_px)
    fig_w = layout.figure_width_px;
    horizontal_margins = layout.left_margin_px + layout.axes_legend_gap_px ...
        + layout.right_margin_px + layout.legend_print_safety_px;
    available_width_px = fig_w - horizontal_margins;
    if isfield(layout, 'axes_width_fraction') && ~isempty(layout.axes_width_fraction)
        layout.axes_width_px = round(available_width_px * layout.axes_width_fraction);
        layout.legend_max_width_px = available_width_px - layout.axes_width_px;
        return
    end
else
    horizontal_span = layout.left_margin_px + layout.axes_width_px + layout.axes_legend_gap_px ...
        + layout.right_margin_px + layout.legend_print_safety_px;
    if isfield(layout, 'legend_max_width_px') && ~isempty(layout.legend_max_width_px)
        legend_budget_px = layout.legend_max_width_px;
    else
        legend_budget_px = layout.legend_min_width_px;
    end
    fig_w = horizontal_span + legend_budget_px;
    layout.figure_width_px = fig_w;
end

horizontal_span = layout.left_margin_px + layout.axes_width_px + layout.axes_legend_gap_px ...
    + layout.right_margin_px + layout.legend_print_safety_px;
layout.legend_max_width_px = fig_w - horizontal_span;
end
