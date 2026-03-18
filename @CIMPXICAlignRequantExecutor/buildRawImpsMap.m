function raw_imps_map = buildRawImpsMap(obj, rawIdentManager)
% Build per-IMP info maps for all raws in a peptide block.
% Input:
%   obj (CIMPXICAlignRequantExecutor)
%       Executor instance
%   rawIdentManager (CIMPRawIdentManager)
%       Per-raw identification store manager
% Output:
%   raw_imps_map (containers.Map)
%       Map: raw_name -> containers.Map(imp_key -> imp_info struct)
%       imp_info fields:
%         - impName (char/string)
%         - impMass (double)
%         - ratio (double vector)
%         - rts (double vector)
%         - intensity (double vector)
%         - charge (double/int)
%         - lowMzBound (double)
%         - highMzBound (double)

[raw_names, raw_ident_stores] = rawIdentManager.getEntries();

groupAggregator = CIMPGroupAggregator(obj.m_ms1_tolerance);
base_groups = CIMPGroup.empty(0, 1);
base_groups = groupAggregator.aggregate(raw_names, raw_ident_stores, [], ...
    @(state, group) append_group(state, group), base_groups);

raw_imps_map = obj.buildRawImpsMapFromBaseGroups(base_groups);
end


function state = append_group(state, group)
group.impRtRanges = [];
if isempty(state)
    state = group;
else
    state(end + 1, 1) = group; %#ok<AGROW>
end
end
