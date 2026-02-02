function runGather(obj)
% Main entry point for summarizing the quantification of various IMPs of a peptide
% Output to a file
% Below each peptide, there are several lines starting with '*', representing the overall information of IMP.
% Each line starting with '*' (IMP line) is followed by several retention time lines ('@' starting lines).
%
% Input:
%   obj (CIMPGatherQuant)
%       Quantification aggregator instance

imp_records = repmat(struct('imp_name','', 'charge',0, 'raw_name','', ...
    'mass_center',0, 'low_mz_bound',0, 'high_mz_bound',0, 'area',0, ...
    'rt_peaks',struct('start',{},'end',{},'ratio',{},'check_label',{})), 0, 1);

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
            % group quant
            [has_nonzero_imp, imp_idx_nonzero, area_imp_final, rt_bound, idx_selected, ratio_each_XIC_peak] = ...
                obj.quant_each_group(keys_raw{idx_keys},...
                group_ratio(charge_group_idxs{idx_ch},:),...
                group_rts(charge_group_idxs{idx_ch},:),...
                group_inten(charge_group_idxs{idx_ch},:),...
                low_mz_bound(idx_ch),high_mz_bound(idx_ch),selected_charge(idx_ch));

            % Only record the non-zero result
            if ~has_nonzero_imp
                % TODO?
                % Show that there is no quantification for this group
                continue;
            end
            for i_iso = 1:length(imp_idx_nonzero)
                rt_peaks = repmat(struct('start',0,'end',0,'ratio',0,'check_label',0), 0, 1);
                for i_peak = 1:size(rt_bound,2)
                    % If the ratio is 0, skip
                    if ratio_each_XIC_peak(i_iso,i_peak) == 0
                        continue;
                    end
                    if i_peak == idx_selected(i_iso)
                        check_label = 1;
                    else
                        check_label = 0;
                    end
                    rt_peaks(end+1,1) = struct(...
                        'start', rt_bound(i_iso,i_peak).start, ...
                        'end', rt_bound(i_iso,i_peak).end, ...
                        'ratio', ratio_each_XIC_peak(i_iso,i_peak), ...
                        'check_label', check_label); %#ok<AGROW>
                end
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
