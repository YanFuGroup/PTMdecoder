function raw_imps_map = buildRawImpsMapFromBaseGroups(~, base_groups)
% Build per-IMP info maps from prebuilt grouped contexts.
% Input:
%   base_groups (CIMPGroup array)
%       grouped contexts for one peptide
% Output:
%   raw_imps_map (containers.Map)
%       Map: raw_name -> containers.Map(imp_key -> imp_info struct)

raw_imps_map = containers.Map('KeyType', 'char', 'ValueType', 'any');
if isempty(base_groups)
    return;
end

for idx_group = 1:numel(base_groups)
    group = base_groups(idx_group);

    if ~isKey(raw_imps_map, group.rawName)
        raw_imps_map(group.rawName) = containers.Map('KeyType', 'char', 'ValueType', 'any');
    end
    raw_map = raw_imps_map(group.rawName);

    for idx_imp = 1:numel(group.impNames)
        key = sprintf('%s|%d', group.impNames{idx_imp}, group.selectedCharge);
        base_key = CIMPQuantRecord.build_id(group.impNames{idx_imp}, group.selectedCharge, group.rawName);
        imp_info = struct( ...
            'impName', group.impNames{idx_imp}, ...
            'impMass', group.impMass(idx_imp), ...
            'ratio', group.ratio(group.chargeGroupIdxs, idx_imp), ...
            'rts', group.rts(group.chargeGroupIdxs, :), ...
            'intensity', group.intensity(group.chargeGroupIdxs, :), ...
            'charge', group.selectedCharge, ...
            'baseKey', base_key, ...
            'lowMzBound', group.lowMzBound, ...
            'highMzBound', group.highMzBound);
        raw_map(key) = imp_info;
    end

    raw_imps_map(group.rawName) = raw_map;
end
end
