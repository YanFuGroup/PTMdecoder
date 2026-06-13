function cfg = fromParamMap(task_param_map)
% Build peptide-level quant service config struct from parameter map.
% Input:
%   task_param_map (containers.Map)
%       parameter key-value map for peptide-level quant service
% Output:
%   cfg (struct)
%       normalized service config
if ~isa(task_param_map, 'containers.Map')
    error('CPeptideQuantServiceConfig:InvalidParamMap', ...
        'Expected task_param_map to be a containers.Map.');
end

cfg = struct();
cfg.mod_file_path = CParamMapUtils.getRequired(task_param_map, CPTMdecoderWorkflowParamKeys.PARAM_MOD_FILE_PATH, 'modification file path', 'CPeptideQuantServiceConfig');
cfg.fixed_mod = CParamMapUtils.getRequired(task_param_map, CPTMdecoderWorkflowParamKeys.PARAM_FIXED_MOD, 'fixed modifications', 'CPeptideQuantServiceConfig');
cfg.variable_mod = CParamMapUtils.getRequired(task_param_map, CPTMdecoderWorkflowParamKeys.PARAM_VARIABLE_MOD, 'variable modifications', 'CPeptideQuantServiceConfig');
cfg.spec_dir_path = CParamMapUtils.getRequired(task_param_map, CPTMdecoderWorkflowParamKeys.PARAM_SPEC_DIR_PATH, 'spectrum directory path', 'CPeptideQuantServiceConfig');

ms1_tol_value = CParamMapUtils.getRequiredNumber(task_param_map, CPTMdecoderWorkflowParamKeys.PARAM_MS1_TOLERANCE_VALUE, 'MS1 tolerance value', 'CPeptideQuantServiceConfig');
ms1_tol_type = strtrim(CParamMapUtils.getRequired(task_param_map, CPTMdecoderWorkflowParamKeys.PARAM_MS1_TOLERANCE_TYPE, 'MS1 tolerance type', 'CPeptideQuantServiceConfig'));
cfg.ms1_tolerance = struct('value', ms1_tol_value, 'isppm', strcmpi(ms1_tol_type, 'PPM'));

cfg.alpha = CParamMapUtils.getRequiredNumber(task_param_map, CPTMdecoderWorkflowParamKeys.PARAM_ALPHA, 'noise filtering threshold', 'CPeptideQuantServiceConfig');
cfg.fasta_file_path = CParamMapUtils.getRequired(task_param_map, CPTMdecoderWorkflowParamKeys.PARAM_FASTA_FILE_PATH, 'FASTA file path', 'CPeptideQuantServiceConfig');
cfg.regular_express = CParamMapUtils.getRequired(task_param_map, CPTMdecoderWorkflowParamKeys.PARAM_REGULAR_EXPRESS, 'regular expression', 'CPeptideQuantServiceConfig');
cfg.filtered_res_file_path = CParamMapUtils.getOptional(task_param_map, CPTMdecoderWorkflowParamKeys.PARAM_FILTERED_RES_FILE_PATH, '', 'CPeptideQuantServiceConfig');
cfg.result_filter_threshold = CParamMapUtils.getRequiredNumber(task_param_map, CPTMdecoderWorkflowParamKeys.PARAM_RESULT_FILTER_THRESHOLD, 'result filter threshold', 'CPeptideQuantServiceConfig');
cfg.output_dir_path = CParamMapUtils.getRequired(task_param_map, CPTMdecoderWorkflowParamKeys.PARAM_OUTPUT_DIR_PATH, 'output directory path', 'CPeptideQuantServiceConfig');

cfg.min_MSMS_num = CParamMapUtils.getOptionalNumber(task_param_map, CPTMdecoderWorkflowParamKeys.PARAM_MIN_MSMS_NUM, 1, 'CPeptideQuantServiceConfig');
cfg.min_xic_nonzero_points = CParamMapUtils.getOptionalNumber(task_param_map, CPTMdecoderWorkflowParamKeys.PARAM_MIN_XIC_NONZERO_POINTS, 5, 'CPeptideQuantServiceConfig');
validateMinXicNonzeroPoints(cfg.min_xic_nonzero_points, 'CPeptideQuantServiceConfig');
cfg.msms_res_path = CParamMapUtils.getOptional(task_param_map, CPTMdecoderWorkflowParamKeys.PARAM_MSMS_RES_PATH, [], 'CPeptideQuantServiceConfig');
cfg.msms_stability_filter = CMsmsStabilityFilterConfig.fromParamMap(task_param_map, 'CPeptideQuantServiceConfig');
end

function validateMinXicNonzeroPoints(value, context)
if ~isscalar(value) || ~isfinite(value) || value < 1 || floor(value) ~= value
    CLogger.error('[%s:InvalidMinXicNonzeroPoints] min_xic_nonzero_points must be a positive integer.', context);
end
end
