function rt_ranges = get_rt_ranges(obj, mod_peptide)
% Get retention time ranges for a modified peptide
% Input:
%   obj (CPepResReader)
%       Peptide result reader instance
%   mod_peptide (1 x 1 char/string)
%       Modified peptide string
% Output:
%   rt_ranges (1 x R struct)
%       fields: rt_start (minutes), rt_end (minutes), ratio, check_label

if isKey(obj.m_pep_res, mod_peptide)
    rt_ranges = obj.m_pep_res(mod_peptide).rt_ranges;
else
    rt_ranges = [];
end
end