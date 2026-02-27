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

cfg = struct();
cfg.mod_file_path = CParamMapUtils.getRequired(task_param_map, CPTMdecoderWorkflowParamKeys.PARAM_MOD_FILE_PATH, 'modification file path', 'CPeptideAlignRequantServiceConfig');
cfg.fixed_mod = CParamMapUtils.getRequired(task_param_map, CPTMdecoderWorkflowParamKeys.PARAM_FIXED_MOD, 'fixed modifications', 'CPeptideAlignRequantServiceConfig');
cfg.variable_mod = CParamMapUtils.getRequired(task_param_map, CPTMdecoderWorkflowParamKeys.PARAM_VARIABLE_MOD, 'variable modifications', 'CPeptideAlignRequantServiceConfig');
cfg.spec_dir_path = CParamMapUtils.getRequired(task_param_map, CPTMdecoderWorkflowParamKeys.PARAM_SPEC_DIR_PATH, 'spectrum directory path', 'CPeptideAlignRequantServiceConfig');

ms1_tol_value = CParamMapUtils.getRequiredNumber(task_param_map, CPTMdecoderWorkflowParamKeys.PARAM_MS1_TOLERANCE_VALUE, ...
    'MS1 tolerance value', 'CPeptideAlignRequantServiceConfig');
ms1_tol_type = strtrim(CParamMapUtils.getRequired(task_param_map, CPTMdecoderWorkflowParamKeys.PARAM_MS1_TOLERANCE_TYPE, ...
    'MS1 tolerance type', 'CPeptideAlignRequantServiceConfig'));
cfg.ms1_tolerance = struct('value', ms1_tol_value, 'isppm', strcmpi(ms1_tol_type, 'PPM'));

cfg.alpha = CParamMapUtils.getRequired(task_param_map, CPTMdecoderWorkflowParamKeys.PARAM_ALPHA, 'noise filtering threshold', 'CPeptideAlignRequantServiceConfig');
cfg.result_filter_threshold = CParamMapUtils.getRequired(task_param_map, CPTMdecoderWorkflowParamKeys.PARAM_RESULT_FILTER_THRESHOLD, 'result filter threshold', 'CPeptideAlignRequantServiceConfig');
cfg.fasta_file_path = CParamMapUtils.getRequired(task_param_map, CPTMdecoderWorkflowParamKeys.PARAM_FASTA_FILE_PATH, 'FASTA file path', 'CPeptideAlignRequantServiceConfig');
cfg.regular_express = CParamMapUtils.getRequired(task_param_map, CPTMdecoderWorkflowParamKeys.PARAM_REGULAR_EXPRESS, 'regular expression', 'CPeptideAlignRequantServiceConfig');
cfg.output_dir_path = CParamMapUtils.getRequired(task_param_map, CPTMdecoderWorkflowParamKeys.PARAM_OUTPUT_DIR_PATH, 'output directory path', 'CPeptideAlignRequantServiceConfig');

cfg.filtered_res_file_path = CParamMapUtils.getOptional(task_param_map, CPTMdecoderWorkflowParamKeys.PARAM_FILTERED_RES_FILE_PATH, '');
cfg.msms_res_path = CParamMapUtils.getOptional(task_param_map, CPTMdecoderWorkflowParamKeys.PARAM_MSMS_RES_PATH, []);
cfg.min_MSMS_num = CParamMapUtils.getOptional(task_param_map, CPTMdecoderWorkflowParamKeys.PARAM_MIN_MSMS_NUM, 1);
end
