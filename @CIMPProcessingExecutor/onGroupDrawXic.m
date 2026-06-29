function state = onGroupDrawXic(obj, state, group)
% Draw XICs for a single group if RT ranges are available
% Input:
%   obj (CIMPProcessingExecutor)
%       processing executor instance
%   state (struct)
%       draw state with fields: dir_save, color_map, legend_map, xic_layout
%   group (CIMPGroup)
%       group payload for one raw and charge
% Output:
%   state (struct)
%       updated draw state

if isempty(group.impRtRanges) || all(cellfun(@isempty, group.impRtRanges))
    return;
end

CIMPQuantPlotter.drawGroup(obj.m_ms12DatasetIO, obj.m_minMSMSnum, group.rawName, ...
    group.ratio(group.chargeGroupIdxs,:), ...
    group.rts(group.chargeGroupIdxs,:), ...
    group.intensity(group.chargeGroupIdxs,:), ...
    group.lowMzBound, group.highMzBound, ...
    group.selectedCharge, group.impRtRanges, ...
    group.impNames, state.dir_save, state.color_map, state.legend_map, ...
    state.xic_layout);
end
