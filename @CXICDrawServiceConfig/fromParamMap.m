function cfg = fromParamMap(task_param_map)
% Build draw-XIC config struct from parameter map.
% Input:
%   task_param_map (containers.Map)
%       parameter key-value map for draw-XIC service
% Output:
%   cfg (struct)
%       normalized service config
if ~isa(task_param_map, 'containers.Map')
    error('CXICDrawServiceConfig:InvalidParamMap', ...
        'Expected task_param_map to be a containers.Map.');
end

cfg = struct();
cfg.mod_file_path = CParamMapUtils.getRequired(task_param_map, CPTMdecoderWorkflowParamKeys.PARAM_MOD_FILE_PATH, 'modification file path', 'CXICDrawServiceConfig');
cfg.fixed_mod = CParamMapUtils.getRequired(task_param_map, CPTMdecoderWorkflowParamKeys.PARAM_FIXED_MOD, 'fixed modifications', 'CXICDrawServiceConfig');
cfg.variable_mod = CParamMapUtils.getRequired(task_param_map, CPTMdecoderWorkflowParamKeys.PARAM_VARIABLE_MOD, 'variable modifications', 'CXICDrawServiceConfig');
cfg.spec_dir_path = CParamMapUtils.getRequired(task_param_map, CPTMdecoderWorkflowParamKeys.PARAM_SPEC_DIR_PATH, 'spectrum directory path', 'CXICDrawServiceConfig');

ms1_tol_value = CParamMapUtils.getRequiredNumber(task_param_map, CPTMdecoderWorkflowParamKeys.PARAM_MS1_TOLERANCE_VALUE, ...
    'MS1 tolerance value', 'CXICDrawServiceConfig');
ms1_tol_type = strtrim(CParamMapUtils.getRequired(task_param_map, CPTMdecoderWorkflowParamKeys.PARAM_MS1_TOLERANCE_TYPE, ...
    'MS1 tolerance type', 'CXICDrawServiceConfig'));
cfg.ms1_tolerance = struct('value', ms1_tol_value, 'isppm', strcmpi(ms1_tol_type, 'PPM'));

cfg.alpha = CParamMapUtils.getRequired(task_param_map, CPTMdecoderWorkflowParamKeys.PARAM_ALPHA, 'noise filtering threshold', 'CXICDrawServiceConfig');
cfg.result_filter_threshold = CParamMapUtils.getRequired(task_param_map, CPTMdecoderWorkflowParamKeys.PARAM_RESULT_FILTER_THRESHOLD, 'result filter threshold', 'CXICDrawServiceConfig');
cfg.output_dir_path = CParamMapUtils.getRequired(task_param_map, CPTMdecoderWorkflowParamKeys.PARAM_OUTPUT_DIR_PATH, 'output directory path', 'CXICDrawServiceConfig');

cfg.checked_peptides_res_path = CParamMapUtils.getOptional(task_param_map, CPTMdecoderWorkflowParamKeys.PARAM_CHECKED_PEPTIDES_RES_PATH, []);
cfg.msms_res_path = CParamMapUtils.getOptional(task_param_map, CPTMdecoderWorkflowParamKeys.PARAM_MSMS_RES_PATH, []);
cfg.min_MSMS_num = CParamMapUtils.getOptional(task_param_map, CPTMdecoderWorkflowParamKeys.PARAM_MIN_MSMS_NUM, 1);
end
