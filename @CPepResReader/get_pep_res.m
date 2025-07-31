function pep_res = get_pep_res(obj)
% Get the peptide result structure
% Output:
%   pep_res - A structure containing the peptide results with fields:
%       peptidoform_name: The modified peptide name
%       charge: The charge state of the peptide
%       dataset_name: The dataset name
%       mean_mz: The mean m/z value
%       lb_mz: The lower bound m/z value
%       ub_mz: The upper bound m/z value
%       quant_value: The quantification value of the peptide
%       rt_ranges: The retention time ranges for the peptide

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