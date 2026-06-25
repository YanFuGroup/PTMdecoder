function layout = getXicLegendLayoutConfig()
% Return configurable geometry for XIC legend layout.
layout = struct();

% Common figure tuning parameters. Adjust these first when preparing final
% figures for publication.
layout.dpi = 300;
layout.figure_width_px = 2000;
layout.figure_height_px = 800;
layout.axes_width_fraction = 0.7;
layout.left_margin_px = 160;
layout.bottom_margin_px = 160;
layout.top_margin_px = 80;
layout.axes_legend_gap_px = 40;
layout.right_margin_px = 30;
layout.legend_chars_per_px = 0.15;
layout.all_font_size = 7;
layout.all_line_width = 1;

% Advanced legend layout safeguards. These usually do not need manual
% changes unless export clipping or extreme label wrapping reappears.
layout.legend_min_width_px = 260;
layout.legend_padding_px = 20;
layout.legend_print_safety_px = 32;
layout.legend_print_width_scale = 1.15;
layout.legend_min_line_chars = 4;
layout.legend_default_icon_width_px = 30;

% Derived horizontal geometry. Keep axes_width_px and legend_max_width_px
% coupled through axes_width_fraction so the fixed figure width is preserved.
available_width_px = layout.figure_width_px ...
    - layout.left_margin_px - layout.axes_legend_gap_px ...
    - layout.right_margin_px - layout.legend_print_safety_px;
layout.axes_width_px = round(available_width_px * layout.axes_width_fraction);
layout.legend_max_width_px = available_width_px - layout.axes_width_px;
end
