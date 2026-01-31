function pep_rtrange_map = get_pep_rtrange_map(obj)
% Get the map of modified peptide to retention time ranges
% Input:
%   obj (CPepResReader)
%       Peptide result reader instance
% Output:
%   pep_rtrange_map (containers.Map)
%       key: modified peptide string; value: rt_ranges struct array
pep_rtrange_map = containers.Map();
if ~isempty(obj.m_pep_res)
    keys = obj.m_pep_res.keys();
    for i = 1:length(keys)
        key = keys{i};
        pep_rtrange_map(key) = obj.m_pep_res(key).rt_ranges;
    end
end
end