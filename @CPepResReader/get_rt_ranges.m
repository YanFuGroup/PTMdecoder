function rt_ranges = get_rt_ranges(obj, mod_peptide)
% Get retention time ranges for a modified peptide
% Input:
%   mod_peptide - The modified peptide string
% Output:
%   rt_ranges - The retention time ranges, ratio and check label for the specified modified peptide
%      rt_start: The start time of the retention time range
%      rt_end: The end time of the retention time range
%      rate: The rate of the retention time range
%      check_label: The check label of the retention time range

if isKey(obj.m_pep_res, mod_peptide)
    rt_ranges = obj.m_pep_res(mod_peptide).rt_ranges;
else
    rt_ranges = [];
end
end