function current_imp_rt_range = buildImpRtRanges(~, group_imp_name, selected_charge, raw_name, pep_rtrange_map)
    % Build retention time ranges for IMPs in a group
    % input:
    %   obj (CIMPGroupAggregator)
    %       the aggregator instance (ignored)
    %   group_imp_name (cell array of strings)
    %       names of IMPs in the current group
    %   selected_charge (1 x 1 double/int)
    %       selected precursor charge state
    %   raw_name (1 x 1 char/string)
    %       name of the raw file
    %   pep_rtrange_map (containers.Map)
    %       map of [modified peptide _+ charge _ raw file name] -> [rt_start, rt_end, check_label]
    % output:
    %   current_imp_rt_range (cell array)
    %       RT ranges for each IMP (each cell: [] or [rt_start, rt_end, check_label])
    if isempty(pep_rtrange_map)
        current_imp_rt_range = [];
        return;
    end

    current_imp_rt_range = cell(length(group_imp_name), 1);
    for idx_imp = 1:length(group_imp_name)
        generated_key = CIMPQuantRecord.build_id(group_imp_name{idx_imp}, selected_charge, raw_name);
        if pep_rtrange_map.isKey(generated_key)
            current_imp_rt_range{idx_imp} = pep_rtrange_map(generated_key);
        end
    end
end