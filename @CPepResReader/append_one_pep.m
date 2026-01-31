function obj = append_one_pep(obj, ...
    key_mod_pep, peptidoform_name, charge, dataset_name, mean_mz, lb_mz, ub_mz, quant_value, rt_ranges, ori_pep_line)
% Add one peptide record
% Input:
%   obj (CPepResReader)
%       Peptide result reader instance
%   key_mod_pep (1 x 1 char/string)
%       modified peptide key (peptide_charge_dataset)
%   peptidoform_name (1 x 1 char/string)
%       peptidoform name
%   charge (1 x 1 double/int)
%       charge state of the peptide
%   dataset_name (1 x 1 char/string)
%       dataset name
%   mean_mz (1 x 1 double) m/z
%       mean m/z value
%   lb_mz (1 x 1 double) m/z
%       lower bound m/z value
%   ub_mz (1 x 1 double) m/z
%       upper bound m/z value
%   quant_value (1 x 1 double)
%       quantification value of the peptide
%   rt_ranges (1 x R struct)
%       retention time ranges (fields: rt_start, rt_end, ratio, check_label)
%   ori_pep_line (1 x 1 char/string)
%       original peptide line from the report file

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