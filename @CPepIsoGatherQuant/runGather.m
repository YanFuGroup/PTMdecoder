function runGather(obj)
% Main entry point for summarizing the quantification of various IMPs of a peptide
% Output to a file
% Below each peptide, there are several lines starting with '*', representing the overall information of IMP.
% Each line starting with '*' (IMP line) is followed by several retention time lines ('@' starting lines).
%
% Input:
%   obj (CPepIsoGatherQuant)
%       Quantification aggregator instance

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
            % group quant
            [has_nonzero_imp, imp_idx_nonzero, area_imp_final, rt_bound, idx_selected, ratio_each_XIC_peak] = ...
                obj.quant_each_group(keys_raw{idx_keys},...
                group_ratio(charge_group_idxs{idx_ch},:),...
                group_rts(charge_group_idxs{idx_ch},:),...
                group_inten(charge_group_idxs{idx_ch},:),...
                low_mz_bound(idx_ch),high_mz_bound(idx_ch),selected_charge(idx_ch));

            % Save to file
            % only write the non zero result
            if ~has_nonzero_imp
                % TODO?
                % Show that there is no quantification for this group
                continue;
            end
            if ~has_written_protein
                has_written_protein = true;
                write_protein_start_position_line(fout, obj.m_prot_names_pos);
            end
            for i_iso = 1:length(imp_idx_nonzero)
                % Write the line of peptide, charge, mass, area
                fprintf(fout, '*\t%s\t+%d\t%s\t%.4f\t%f\t%f\t%f\n', ...
                    group_imp_name{imp_idx_nonzero(i_iso)}, ...
                    selected_charge(idx_ch),...
                    keys_raw{idx_keys},...
                    mean([low_mz_bound(idx_ch),high_mz_bound(idx_ch)]),...
                    low_mz_bound(idx_ch),...
                    high_mz_bound(idx_ch),...
                    area_imp_final(i_iso,1));
                % Write the line of rt range, ratio, check label
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
                    fprintf(fout, '@\t%f\t%f\t%f\t%d\n',...
                        rt_bound(i_iso,i_peak).start,...
                        rt_bound(i_iso,i_peak).end, ...
                        ratio_each_XIC_peak(i_iso,i_peak), ...
                        check_label);
                end
            end
        end
    end
end
fclose(fout);

end
