function labels = prepareXicLegendLabels(labels, layout, max_line_chars)
% Prepare legend labels without TeX semantics or redundant prefixes.
if nargin < 3 || isempty(max_line_chars)
    max_line_chars = get_legend_max_line_chars(layout);
end
for idx = 1:numel(labels)
    label = char(string(labels{idx}));
    label = regexprep(label, '^\s*XIC\s+of\s+', '', 'ignorecase');
    if numel(label) > max_line_chars
        label = wrap_legend_label(label, max_line_chars);
    end
    labels{idx} = label;
end
end

function max_line_chars = get_legend_max_line_chars(layout)
if isfield(layout, 'legend_default_icon_width_px') && ~isempty(layout.legend_default_icon_width_px)
    icon_width_px = layout.legend_default_icon_width_px;
else
    icon_width_px = 30;
end
text_width_px = layout.legend_max_width_px - layout.legend_padding_px - icon_width_px;
text_width_px = max(1, text_width_px);
if isfield(layout, 'legend_chars_per_px') && ~isempty(layout.legend_chars_per_px)
    chars_per_px = layout.legend_chars_per_px;
else
    chars_per_px = 0.13;
end
max_line_chars = floor(text_width_px * chars_per_px);
if isfield(layout, 'legend_min_line_chars') && ~isempty(layout.legend_min_line_chars)
    max_line_chars = max(layout.legend_min_line_chars, max_line_chars);
else
    max_line_chars = max(1, max_line_chars);
end
end

function wrapped = wrap_legend_label(label, max_line_chars)
% Wrap a legend label by fixed character count while preserving all characters.
max_line_chars = max(1, max_line_chars);
lines = {};
for start_idx = 1:max_line_chars:numel(label)
    end_idx = min(start_idx + max_line_chars - 1, numel(label));
    lines{end + 1} = label(start_idx:end_idx); %#ok<AGROW>
end
wrapped = strjoin(lines, newline);
end
