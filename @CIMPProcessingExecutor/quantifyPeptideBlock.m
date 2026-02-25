function block = quantifyPeptideBlock(obj, prot_names_pos, rawIdentManager)
% Build a quantification block for one peptide
% Input:
%   obj (CIMPProcessingExecutor)
%       processing executor instance
%   prot_names_pos (P x 2 cell)
%       protein name and start position pairs
%   rawIdentManager (CIMPRawIdentManager)
%       per-raw identification store manager
% Output:
%   block (CIMPQuantBlock or empty)
%       protein block with IMP records, empty if no records

imp_records = CIMPQuantRecord.empty(0,1);
[raw_names, raw_ident_stores] = rawIdentManager.getEntries();

groupAggregator = CIMPGroupAggregator(obj.m_ms1_tolerance);
imp_records = groupAggregator.aggregate(raw_names, raw_ident_stores, [], ...
    @(state, group) obj.onGroupQuant(state, group), imp_records);

if isempty(imp_records)
    block = CIMPQuantBlock.empty(0,1);
else
    block = CIMPQuantBlock(prot_names_pos, imp_records);
end
end
