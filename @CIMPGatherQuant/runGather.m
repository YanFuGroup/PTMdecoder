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


% Do the same operation for gathered PSM for every raw file
imp_records = iterate_imp_groups(obj, [], @handle_group_charge, imp_records);
CIMPGatherWriter.write_imp_group_block(obj.m_outputPath, obj.m_prot_names_pos, imp_records);
end



function imp_records = handle_group_charge(imp_records, obj, raw_name, group_imp_name, group_ratio, group_rts, group_inten, ...
        ~, low_mz_bound, high_mz_bound, selected_charge, charge_group_idxs, ~)
% group quant
[has_nonzero_imp, imp_idx_nonzero, area_imp_final, rt_bound, idx_selected, ratio_each_XIC_peak] = ...
    obj.quant_each_group(raw_name,...
    group_ratio(charge_group_idxs,:),...
    group_rts(charge_group_idxs,:),...
    group_inten(charge_group_idxs,:),...
    low_mz_bound,high_mz_bound,selected_charge);

% Only record the non-zero result
if ~has_nonzero_imp
    % TODO?
    % Show that there is no quantification for this group
    return;
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
        'charge', selected_charge, ...
        'raw_name', raw_name, ...
        'mass_center', mean([low_mz_bound,high_mz_bound]), ...
        'low_mz_bound', low_mz_bound, ...
        'high_mz_bound', high_mz_bound, ...
        'area', area_imp_final(i_iso,1), ...
        'rt_peaks', rt_peaks); %#ok<AGROW>
end
end
