function rerunGather_quant(obj,pep_rtrange_map)
% Re-quantification for gathered peptides using manually-checked rt range
% Input:
%   obj (CPepIsoGatherQuant)
%       Quantification aggregator instance
%   pep_rtrange_map (containers.Map)
%       map of [modified peptide _ charge _ raw file name] -> [rt_start, rt_end, check_label]

fout = fopen(obj.m_outputPath,'a');
if fout == -1
    error(['Cannot open the the report file ',obj.m_outputPath]);
end
has_written_protein = false;

% Delete useless raws which are empty or blank.
obj = obj.delUselessRaws();

% Do the same operation for gathered PSM for every raw file
keys_raw = obj.m_mapRawNames.keys;
for idx_keys = 1:obj.m_mapRawNames.Count
    idx_r = obj.m_mapRawNames(keys_raw{idx_keys});

    % Cluster the IMPs according to their masses
    group_idxs = cluster_imps_by_mass(obj.m_IMPMass{idx_r},obj.m_ms1_tolerance);

    % Quantify the IMPs in each group
    for idx_g = 1:length(group_idxs)
        group_imp_name = obj.m_cstrIMPNames{idx_r}(group_idxs{idx_g});
        group_ratio = obj.m_ratioMatrix{idx_r}(:,group_idxs{idx_g});
        idxs_rt_inten = find(sum(group_ratio,2));
        group_ratio = group_ratio(idxs_rt_inten,:);
        group_rts = obj.m_curRts{idx_r}(idxs_rt_inten);
        group_inten = obj.m_curIntens{idx_r}(idxs_rt_inten);
        group_imp_mass = obj.m_IMPMass{idx_r}(group_idxs{idx_g});
        group_charge = obj.m_curCharge{idx_r}(idxs_rt_inten);
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

            % Save to file
            % only write the non zero result
            if ~has_nonzero_imp
                continue;
            end
            if ~has_written_protein
                has_written_protein = true;
                write_protein_start_position_line(fout, obj.m_prot_names_pos);
            end
            for i_iso = 1:length(imp_idx_nonzero)
                fprintf(fout, '*\t%s\t+%d\t%s\t%.4f\t%f\t%f\t%f\n', ...
                    group_imp_name{imp_idx_nonzero(i_iso)}, ...
                    selected_charge(idx_ch),...
                    keys_raw{idx_keys},...
                    mean([low_mz_bound(idx_ch),high_mz_bound(idx_ch)]),...
                    low_mz_bound(idx_ch),...
                    high_mz_bound(idx_ch),...
                    area_imp_final(i_iso,1));
                fprintf(fout, '@\t%f\t%f\t%f\t%d\n', ...
                    rt_bound(i_iso).start, ...
                    rt_bound(i_iso).end, ...
                    ratio_each_XIC_peak(i_iso), ...
                    max_label(i_iso));
            end
        end
    end
end
fclose(fout);
end