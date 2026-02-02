function state = iterate_imp_groups(obj, pep_rtrange_map, on_group_charge, state)
% Iterate over raw->group->charge and invoke a callback with state accumulation
% Input:
%   obj (CIMPGatherQuant)
%   pep_rtrange_map (containers.Map or [])
%       map of [modified peptide _ charge _ raw file name] -> [rt_start, rt_end, check_label]
%   on_group_charge (function_handle)
%       callback for each (raw, group, charge), signature:
%       state = on_group_charge(state, obj, group)
%   state (any)
%       accumulator passed through each callback

if nargin < 2
    pep_rtrange_map = [];
end
if nargin < 4
    state = [];
end

% Do the same operation for gathered PSM for every raw file
keys_raw = obj.m_mapRawNames.keys;
for idx_keys = 1:obj.m_mapRawNames.Count
    idx_r = obj.m_mapRawNames(keys_raw{idx_keys});
    raw = obj.get_raw(idx_r);

    if raw.length <= 0
        error('CIMPGatherQuant:InvalidRawLength', ...
            'raw.length must be > 0 (raw: %s).', keys_raw{idx_keys});
    end
    ratio_cols = size(raw.ratioMatrix,2);
    if length(raw.impNames) ~= ratio_cols || length(raw.impMass) ~= ratio_cols
        error('CIMPGatherQuant:InvalidRawRatioColumns', ...
            'impNames/impMass length must match ratioMatrix columns (raw: %s).', ...
            keys_raw{idx_keys});
    end

    % Cluster the IMPs according to their masses
    group_idxs = CIMPGrouper.cluster_by_mass(raw.impMass,obj.m_ms1_tolerance);

    % Quantify the IMPs in each group
    for idx_g = 1:length(group_idxs)
        group_imp_name = raw.impNames(group_idxs{idx_g});
        cur_ratio = raw.ratioMatrix(1:raw.length,group_idxs{idx_g});
        cur_rts = raw.curRts(1:raw.length);
        cur_inten = raw.curIntens(1:raw.length);
        cur_charge = raw.curCharge(1:raw.length);
        idxs_rt_inten = find(sum(cur_ratio,2));
        group_ratio = cur_ratio(idxs_rt_inten,:);
        group_rts = cur_rts(idxs_rt_inten);
        group_inten = cur_inten(idxs_rt_inten);
        group_imp_mass = raw.impMass(group_idxs{idx_g});
        group_charge = cur_charge(idxs_rt_inten);
        [low_mz_bound, high_mz_bound, selected_charge, charge_group_idxs] = ...
            obj.get_mz_bound(group_imp_mass,group_charge);

        for idx_ch = 1:length(selected_charge)
            current_imp_rt_range = [];
            if ~isempty(pep_rtrange_map)
                % Get retention time range for each IMP
                current_imp_rt_range = cell(length(group_imp_name),1);
                for idx_imp = 1:length(group_imp_name)
                    generated_key = [group_imp_name{idx_imp},'_+', ...
                        num2str(selected_charge(idx_ch)), '_', keys_raw{idx_keys}];
                    if pep_rtrange_map.isKey(generated_key)
                        current_imp_rt_range{idx_imp} = pep_rtrange_map(generated_key);
                    end
                end
            end

            group = CIMPGroup(keys_raw{idx_keys}, group_imp_name, group_imp_mass, ...
                group_ratio, group_rts, group_inten, group_charge, ...
                low_mz_bound(idx_ch), high_mz_bound(idx_ch), selected_charge(idx_ch), ...
                charge_group_idxs{idx_ch}, current_imp_rt_range);

            state = on_group_charge(state, obj, group);
        end
    end
end

end
