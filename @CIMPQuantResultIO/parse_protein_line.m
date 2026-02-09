function protein_name_pos = parse_protein_line(strline)
% Parse a protein line (not starting with '*' or '@')
% Line format: protein1,start1;protein2,start2;...
% Semicolon-separated pairs, each pair is comma-separated protein name and start position
protein_name_pos = {};
segments = regexp(strline, ';', 'split');
for idx_seg = 1:numel(segments)
    token = strtrim(segments{idx_seg});
    if isempty(token)
        continue;
    end
    pair = regexp(token, ',', 'split');
    if numel(pair) < 2
        continue;
    end
    protein_name = strtrim(pair{1});
    start_pos = str2double(strtrim(pair{2}));
    protein_name_pos(end+1,1:2) = {protein_name, start_pos}; %#ok<AGROW>
end
end
