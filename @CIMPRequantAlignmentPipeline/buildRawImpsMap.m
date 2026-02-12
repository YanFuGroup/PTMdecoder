function raw_imps_map = buildRawImpsMap(obj, rawIdentManager)
% Build per-IMP info maps for all raws in a peptide block.
% Input:
%   obj (CIMPRequantAlignmentPipeline)
%       Pipeline instance
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

raw_imps_map = containers.Map('KeyType', 'char', 'ValueType', 'any');
[raw_names, raw_ident_stores] = rawIdentManager.getEntries();

groupAggregator = CIMPGroupAggregator(obj.m_ms1_tolerance);
groupAggregator.aggregate(raw_names, raw_ident_stores, [], ...
    @(state, group) add_imp_info(state, group), raw_imps_map);
end

function state = add_imp_info(state, group)
if ~isKey(state, group.rawName)
    state(group.rawName) = containers.Map('KeyType', 'char', 'ValueType', 'any');
end
raw_map = state(group.rawName);
for idx_imp = 1:numel(group.impNames)
    key = sprintf('%s|%d', group.impNames{idx_imp}, group.selectedCharge);
    imp_info = struct( ...
        'impName', group.impNames{idx_imp}, ...
        'impMass', group.impMass(idx_imp), ...
        'ratio', group.ratio(group.chargeGroupIdxs, idx_imp), ...
        'rts', group.rts(group.chargeGroupIdxs, :), ...
        'intensity', group.intensity(group.chargeGroupIdxs, :), ...
        'charge', group.selectedCharge, ...
        'lowMzBound', group.lowMzBound, ...
        'highMzBound', group.highMzBound);
    raw_map(key) = imp_info;
end
state(group.rawName) = raw_map;
end
