function block = buildRequantBlock(obj, prot_names_pos, rawManager, pep_rtrange_map)
% Build a re-quantification block using checked RT ranges
% Input:
%   obj (CIMPProcessingPipeline)
%       processing pipeline instance
%   prot_names_pos (P x 2 cell)
%       protein name and start position pairs
%   rawManager (CIMPRawIdentManager)
%       per-raw identification store manager
%   pep_rtrange_map (containers.Map)
%       map of [modified peptide _ charge _ raw file name] -> [rt_start, rt_end, check_label]
% Output:
%   block (CIMPQuantBlock or empty)
%       protein block with IMP records, empty if no records

imp_records = CIMPQuantRecord.empty(0,1);
[raw_names, raw_ident_stores] = rawManager.getEntries();

groupAggregator = CIMPGroupAggregator(obj.m_ms1_tolerance);
imp_records = groupAggregator.aggregate(raw_names, raw_ident_stores, pep_rtrange_map, ...
    @(state, group) obj.appendRequantGroupRecord(state, group), imp_records);

if isempty(imp_records)
    block = CIMPQuantBlock.empty(0,1);
else
    block = CIMPQuantBlock(prot_names_pos, imp_records);
end
end
