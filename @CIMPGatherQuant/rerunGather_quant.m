function rerunGather_quant(obj,pep_rtrange_map)
% Re-quantification for gathered peptides using manually-checked rt range
% Input:
%   obj (CIMPGatherQuant)
%       Quantification aggregator instance
%   pep_rtrange_map (containers.Map)
%       map of [modified peptide _ charge _ raw file name] -> [rt_start, rt_end, check_label]

imp_records = repmat(struct('imp_name','', 'charge',0, 'raw_name','', ...
    'mass_center',0, 'low_mz_bound',0, 'high_mz_bound',0, 'area',0, ...
    'rt_peaks',struct('start',{},'end',{},'ratio',{},'check_label',{})), 0, 1);


% Do the same operation for gathered PSM for every raw file
keys_raw = obj.m_mapRawNames.keys;
for idx_keys = 1:obj.m_mapRawNames.Count
    idx_r = obj.m_mapRawNames(keys_raw{idx_keys});
    raw = obj.get_raw(idx_r);

    % Cluster the IMPs according to their masses
    group_idxs = cluster_imps_by_mass(raw.impMass,obj.m_ms1_tolerance);

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
            % Get retention time range for each IMP
            current_imp_rt_range = cell(length(group_imp_name),1);
            for idx_imp = 1:length(group_imp_name)
                generated_key = [group_imp_name{idx_imp},'_+', ...
                    num2str(selected_charge(idx_ch)), '_', keys_raw{idx_keys}];
                if pep_rtrange_map.isKey(generated_key)
                    current_imp_rt_range{idx_imp} = pep_rtrange_map(generated_key);
                end
            end
            if all(cellfun(@isempty,current_imp_rt_range))
                % All of the IMPs are removed in manual checking
                continue;
            end

            % group quant
            [has_nonzero_imp, imp_idx_nonzero, area_imp_final, rt_bound, max_label, ratio_each_XIC_peak] = ...
                obj.requant_each_group(keys_raw{idx_keys},...
                group_ratio(charge_group_idxs{idx_ch},:),...
                group_rts(charge_group_idxs{idx_ch},:),...
                group_inten(charge_group_idxs{idx_ch},:),...
                low_mz_bound(idx_ch),high_mz_bound(idx_ch), ...
                selected_charge(idx_ch),current_imp_rt_range);

            % Only record the non-zero result
            if ~has_nonzero_imp
                continue;
            end
            for i_iso = 1:length(imp_idx_nonzero)
                rt_peaks = struct(...
                    'start', rt_bound(i_iso).start, ...
                    'end', rt_bound(i_iso).end, ...
                    'ratio', ratio_each_XIC_peak(i_iso), ...
                    'check_label', max_label(i_iso));
                imp_records(end+1,1) = struct(...
                    'imp_name', group_imp_name{imp_idx_nonzero(i_iso)}, ...
                    'charge', selected_charge(idx_ch), ...
                    'raw_name', keys_raw{idx_keys}, ...
                    'mass_center', mean([low_mz_bound(idx_ch),high_mz_bound(idx_ch)]), ...
                    'low_mz_bound', low_mz_bound(idx_ch), ...
                    'high_mz_bound', high_mz_bound(idx_ch), ...
                    'area', area_imp_final(i_iso,1), ...
                    'rt_peaks', rt_peaks); %#ok<AGROW>
            end
        end
    end
end
CIMPGatherWriter.write_imp_group_block(obj.m_outputPath, obj.m_prot_names_pos, imp_records);
end