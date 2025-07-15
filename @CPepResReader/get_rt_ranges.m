function rt_ranges = get_rt_ranges(obj, mod_peptide)
% Get retention time ranges for a modified peptide
% Input:
%   mod_peptide - The modified peptide string
% Output:
%   rt_ranges - The retention time ranges for the specified modified peptide

if isKey(obj.m_pep_res, mod_peptide)
    rt_ranges = obj.m_pep_res(mod_peptide).rt_ranges;
else
    rt_ranges = [];
end
end