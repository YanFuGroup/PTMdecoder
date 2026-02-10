function block = buildQuantBlock(obj, prot_names_pos, rawManager, ms12DatasetIO)
% Build the quantification block for one peptide
% Input:
%   obj (CPepNormalization)
%       Normalization processor instance
%   prot_names_pos (1 x 2 cell)
%       protein name and start position pairs
%   rawManager (CIMPRawIdentManager)
%       per-raw identification store manager
%   ms12DatasetIO (object)
%       MS1/MS2 dataset IO instance
% Output:
%   block (CIMPQuantBlock or empty)
%       Protein block with IMP records, empty if no records

imp_records = CIMPQuantRecord.empty(0,1);
minMSMSnum = 1;

[raw_names, raw_ident_stores] = rawManager.getEntries();
groupAggregator = CIMPGroupAggregator(obj.ms1_tolerance);
imp_records = groupAggregator.aggregate(raw_names, raw_ident_stores, [], ...
    @(state, group) obj.appendQuantGroupRecord(state, group, ms12DatasetIO, minMSMSnum), imp_records);

if isempty(imp_records)
    block = CIMPQuantBlock.empty(0,1);
else
    block = CIMPQuantBlock(prot_names_pos, imp_records);
end
end
