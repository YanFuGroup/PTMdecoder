function pep_res = get_pep_res(obj)
% Get the peptide result structure
% Input:
%   obj (CPepResReader)
%       Peptide result reader instance
% Output:
%   pep_res (1 x N struct)
%       fields: peptidoform_name, charge, dataset_name, mean_mz, lb_mz, ub_mz,
%       quant_value, rt_ranges
%       peptidoform_name: modified peptide name
%       charge: charge state
%       dataset_name: dataset (raw/mgf) name
%       mean_mz: mean m/z value
%       lb_mz: lower bound m/z
%       ub_mz: upper bound m/z
%       quant_value: quantified peak area
%       rt_ranges: struct array of RT ranges (rt_start, rt_end, ratio, check_label)

pep_res = struct('peptidoform_name', {}, 'charge', {}, 'dataset_name', {}, ...
    'mean_mz', {}, 'lb_mz', {}, 'ub_mz', {}, 'quant_value', {}, 'rt_ranges', {});
keys = obj.m_pep_res.keys;
for i = 1:length(keys)
    key = keys{i};
    value = obj.m_pep_res(key);
    
    pep_res(i).peptidoform_name = value.peptidoform_name;
    pep_res(i).charge = value.charge;
    pep_res(i).dataset_name = value.dataset_name;
    pep_res(i).mean_mz = value.mean_mz;
    pep_res(i).lb_mz = value.lb_mz;
    pep_res(i).ub_mz = value.ub_mz;
    pep_res(i).quant_value = value.quant_value;
    pep_res(i).rt_ranges = value.rt_ranges;

end