function obj = append_one_pep(obj, ...
    key_mod_pep, peptidoform_name, charge, dataset_name, mean_mz, lb_mz, ub_mz, quant_value, rt_ranges, ori_pep_line)
% Add one rt ranges to the pep_rtrange struct
% Input:
%   pep_rtrange_map
%       the modified peptide -> retention time range structure
%   key_mod_pep
%       the modified peptide name
%   peptidoform_name
%       the name of the peptidoform
%   charge
%       the charge state of the peptide
%   dataset_name
%       the dataset name
%   mean_mz
%       the mean m/z value
%   lb_mz
%       the lower bound m/z value
%   ub_mz
%       the upper bound m/z value
%   rt_ranges
%       the retention time ranges structure
%   ori_pep_line
%       the original peptide line from the report_peptide_all_checked file

% Skip the first line calling or empty-range peptide.
if isempty(key_mod_pep) || isempty(rt_ranges)
    return;
end

% Create a new peptide record
obj.m_pep_res(key_mod_pep) = struct(...
    'peptidoform_name', peptidoform_name, ...
    'charge', charge, ...
    'dataset_name', dataset_name, ...
    'mean_mz', mean_mz, ...
    'lb_mz', lb_mz, ...
    'ub_mz', ub_mz, ...
    'quant_value', quant_value, ...
    'rt_ranges', rt_ranges, ...
    'ori_pep_line', ori_pep_line ...
);
end