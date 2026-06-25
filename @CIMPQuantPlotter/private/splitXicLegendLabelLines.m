function lines = splitXicLegendLabelLines(label)
% Split a legend label into non-empty display lines.
if contains(label, newline)
    lines = strsplit(label, newline);
else
    lines = {label};
end
lines = lines(~cellfun(@isempty, strtrim(lines)));
if isempty(lines)
    lines = {''};
end
end
