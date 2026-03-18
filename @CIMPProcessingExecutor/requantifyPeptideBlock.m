function block = requantifyPeptideBlock(obj, prot_names_pos, rawIdentManager, pep_rtrange_map, base_groups)
% Build a re-quantification block using checked RT ranges
% Input:
%   obj (CIMPProcessingExecutor)
%       processing executor instance
%   prot_names_pos (P x 2 cell)
%       protein name and start position pairs
%   rawIdentManager (CIMPRawIdentManager)
%       per-raw identification store manager
%   pep_rtrange_map (containers.Map)
%       map of [modified peptide _ charge _ raw file name] -> [rt_start, rt_end, check_label]
%   base_groups (CIMPGroup array, optional)
%       prebuilt grouped contexts for this peptide
% Output:
%   block (CIMPQuantBlock or empty)
%       protein block with IMP records, empty if no records

if nargin < 5 || isempty(base_groups)
    base_groups = obj.buildBaseGroups(rawIdentManager);
end

imp_records = CIMPQuantRecord.empty(0,1);
for idx_group = 1:numel(base_groups)
    group = base_groups(idx_group);
    group.impRtRanges = buildImpRtRangesFromMap(...
        group.impNames, group.selectedCharge, group.rawName, pep_rtrange_map);
    imp_records = obj.onGroupRequant(imp_records, group);
end

if isempty(imp_records)
    block = CIMPQuantBlock.empty(0,1);
else
    block = CIMPQuantBlock(prot_names_pos, imp_records);
end
end


function current_imp_rt_range = buildImpRtRangesFromMap(group_imp_name, selected_charge, raw_name, pep_rtrange_map)
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
