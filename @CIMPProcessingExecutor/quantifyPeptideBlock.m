function block = quantifyPeptideBlock(obj, prot_names_pos, rawIdentManager, base_groups)
% Build a quantification block for one peptide
% Input:
%   obj (CIMPProcessingExecutor)
%       processing executor instance
%   prot_names_pos (P x 2 cell)
%       protein name and start position pairs
%   rawIdentManager (CIMPRawIdentManager)
%       per-raw identification store manager
%   base_groups (CIMPGroup array, optional)
%       prebuilt grouped contexts for this peptide
% Output:
%   block (CIMPQuantBlock or empty)
%       protein block with IMP records, empty if no records

if nargin < 4 || isempty(base_groups)
    base_groups = obj.buildBaseGroups(rawIdentManager);
end

imp_records = CIMPQuantRecord.empty(0,1);
for idx_group = 1:numel(base_groups)
    imp_records = obj.onGroupQuant(imp_records, base_groups(idx_group));
end

if isempty(imp_records)
    block = CIMPQuantBlock.empty(0,1);
else
    block = CIMPQuantBlock(prot_names_pos, imp_records);
end
end
