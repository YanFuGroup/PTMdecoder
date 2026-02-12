function pairs = getAlignmentPairs(obj, raw_names)
% Build reference-run alignment pairs.
% Input:
%   obj (ReferenceRunAlignStrategy)
%       Strategy instance
%   raw_names (1 x N cell)
%       Raw names present in the report
% Output:
%   pairs (M x 2 cell)
%       [ref_raw, target_raw] pairs

if isempty(raw_names)
    pairs = cell(0,2);
    return;
end

if isempty(obj.m_reference_raw)
    obj.m_reference_raw = raw_names{1};
end

pairs = cell(0,2);
for idx = 1:numel(raw_names)
    if ~strcmp(raw_names{idx}, obj.m_reference_raw)
        pairs(end+1, :) = {obj.m_reference_raw, raw_names{idx}};    %#ok<AGROW>
    end
end
end
