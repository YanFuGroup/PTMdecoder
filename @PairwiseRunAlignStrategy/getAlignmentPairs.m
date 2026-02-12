function pairs = getAlignmentPairs(obj, raw_names)
% Validate and return the user-defined alignment pairs.
% Input:
%   obj (PairwiseRunAlignStrategy)
%       Strategy instance with m_pair_list
%   raw_names (1 x N cell)
%       Raw names present in the report
% Output:
%   pairs (M x 2 cell)
%       [ref_raw, target_raw] pairs

pairs = obj.m_pair_list;
if isempty(pairs)
    return;
end
if size(pairs,2) ~= 2
    error('PairwiseRunAlignStrategy:InvalidPairs', ...
        'pair_list must be an N x 2 cell array.');
end

raw_set = containers.Map(raw_names, num2cell(1:numel(raw_names)));
for idx = 1:size(pairs,1)
    if ~isKey(raw_set, pairs{idx,1}) || ~isKey(raw_set, pairs{idx,2})
        error('PairwiseRunAlignStrategy:UnknownRaw', ...
            'Unknown raw in pair: %s -> %s', pairs{idx,1}, pairs{idx,2});
    end
end
end
