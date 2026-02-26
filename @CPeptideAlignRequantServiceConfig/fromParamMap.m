function cfg = fromParamMap(task_param_map)
% Build peptide align-requant config struct from parameter map.
% Input:
%   task_param_map (containers.Map)
%       parameter key-value map for peptide align-requant service
% Output:
%   cfg (struct)
%       normalized service config
if ~isa(task_param_map, 'containers.Map')
    error('CPeptideAlignRequantServiceConfig:InvalidParamMap', ...
        'Expected task_param_map to be a containers.Map.');
end

required_keys = {'mod_file_path', 'fixed_mod', 'variable_mod', ...
    'spec_dir_path', 'ms1_tolerance', 'alpha', 'result_filter_threshold', ...
    'fasta_file_path', 'regular_express', 'output_dir_path'};
for i_key = 1:length(required_keys)
    key_name = required_keys{i_key};
    if ~task_param_map.isKey(key_name)
        error('CPeptideAlignRequantServiceConfig:MissingConfigField', ...
            'Missing required config field: %s', key_name);
    end
end

cfg = struct();
cfg.mod_file_path = task_param_map('mod_file_path');
cfg.fixed_mod = task_param_map('fixed_mod');
cfg.variable_mod = task_param_map('variable_mod');
cfg.spec_dir_path = task_param_map('spec_dir_path');
cfg.ms1_tolerance = task_param_map('ms1_tolerance');
cfg.alpha = task_param_map('alpha');
cfg.result_filter_threshold = task_param_map('result_filter_threshold');
cfg.fasta_file_path = task_param_map('fasta_file_path');
cfg.regular_express = task_param_map('regular_express');
cfg.output_dir_path = task_param_map('output_dir_path');

if task_param_map.isKey('filtered_res_file_path')
    cfg.filtered_res_file_path = task_param_map('filtered_res_file_path');
else
    cfg.filtered_res_file_path = '';
end
if task_param_map.isKey('msms_res_path')
    cfg.msms_res_path = task_param_map('msms_res_path');
else
    cfg.msms_res_path = [];
end
if task_param_map.isKey('min_MSMS_num')
    cfg.min_MSMS_num = task_param_map('min_MSMS_num');
else
    cfg.min_MSMS_num = 1;
end
end
