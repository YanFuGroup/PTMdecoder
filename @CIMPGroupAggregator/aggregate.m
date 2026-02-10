function state = aggregate(obj, raw_names, raw_ident_stores, pep_rtrange_map, handle_group_charge, state)
    % Aggregate IMP groups across raw files and invoke callback
    % input:
    %   obj (CIMPGroupAggregator)
    %       the aggregator instance
    %   raw_names (1 x N cell)
    %       names of the raw files
    %   raw_ident_stores (1 x N cell)
    %       raw identification stores for each raw
    %   pep_rtrange_map (containers.Map or [])
    %       map of [modified peptide _ charge _ raw file name] -> [rt_start, rt_end, check_label]
    %   handle_group_charge (function_handle)
    %       callback function with signature state = f(state, group)
    %   state (any)
    %       accumulator state passed through callbacks
    % output:
    %   state (any)
    %       updated accumulator state after processing all groups
    if nargin < 4
        pep_rtrange_map = [];
    end
    if nargin < 6
        state = [];
    end

    for idx_raw = 1:numel(raw_names)
        raw_name = raw_names{idx_raw};
        raw = raw_ident_stores{idx_raw};

        if raw.length <= 0
            error('CIMPGroupAggregator:InvalidRawLength', ...
                'raw.length must be > 0 (raw: %s).', raw_name);
        end
        ratio_cols = size(raw.ratioMatrix,2);
        if length(raw.impNames) ~= ratio_cols || length(raw.impMass) ~= ratio_cols
            error('CIMPGroupAggregator:InvalidRawRatioColumns', ...
                'impNames/impMass length must match ratioMatrix columns (raw: %s).', ...
                raw_name);
        end

        % Cluster the IMPs according to their masses
        group_idxs = obj.cluster_by_mass(raw.impMass);

        % Aggregate for each group
        for idx_g = 1:length(group_idxs)
            group_imp_name = raw.impNames(group_idxs{idx_g});
            cur_ratio = raw.ratioMatrix(1:raw.length, group_idxs{idx_g});
            cur_rts = raw.curRts(1:raw.length);
            cur_inten = raw.curIntens(1:raw.length);
            cur_charge = raw.curCharge(1:raw.length);
            idxs_rt_inten = find(sum(cur_ratio, 2));
            group_ratio = cur_ratio(idxs_rt_inten, :);
            group_rts = cur_rts(idxs_rt_inten);
            group_inten = cur_inten(idxs_rt_inten);
            group_imp_mass = raw.impMass(group_idxs{idx_g});
            group_charge = cur_charge(idxs_rt_inten);
            [low_mz_bound, high_mz_bound, selected_charge, charge_group_idxs] = ...
                obj.getMzBound(group_imp_mass, group_charge);

            for idx_ch = 1:length(selected_charge)
                current_imp_rt_range = obj.buildImpRtRanges(...
                    group_imp_name, selected_charge(idx_ch), raw_name, pep_rtrange_map);

                group = CIMPGroup(raw_name, group_imp_name, group_imp_mass, ...
                    group_ratio, group_rts, group_inten, group_charge, ...
                    low_mz_bound(idx_ch), high_mz_bound(idx_ch), selected_charge(idx_ch), ...
                    charge_group_idxs{idx_ch}, current_imp_rt_range);

                state = handle_group_charge(state, group);
            end
        end
    end
end