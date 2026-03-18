function base_groups = buildBaseGroups(obj, rawIdentManager)
% Build reusable grouped contexts for one peptide.
% Input:
%   obj (CIMPProcessingExecutor)
%       processing executor instance
%   rawIdentManager (CIMPRawIdentManager)
%       per-raw identification store manager
% Output:
%   base_groups (CIMPGroup array)
%       grouped contexts without requant RT constraints

if isempty(rawIdentManager)
    base_groups = CIMPGroup.empty(0, 1);
    return;
end

[raw_names, raw_ident_stores] = rawIdentManager.getEntries();
groupAggregator = CIMPGroupAggregator(obj.m_ms1_tolerance);
base_groups = CIMPGroup.empty(0, 1);
base_groups = groupAggregator.aggregate(raw_names, raw_ident_stores, [], ...
    @(state, group) append_group(state, group), base_groups);
end


function state = append_group(state, group)
group.impRtRanges = [];
if isempty(state)
    state = group;
else
    state(end + 1, 1) = group; %#ok<AGROW>
end
end
