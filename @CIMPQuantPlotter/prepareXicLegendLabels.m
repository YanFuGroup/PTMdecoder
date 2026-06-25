function labels = prepareXicLegendLabels(labels, layout, max_line_chars)
% Prepare legend labels without TeX semantics or redundant prefixes.
if nargin < 3 || isempty(max_line_chars)
    max_line_chars = getXicLegendMaxLineChars(layout);
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
