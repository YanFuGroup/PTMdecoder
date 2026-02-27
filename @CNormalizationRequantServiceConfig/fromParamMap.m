function cfg = fromParamMap(task_param_map)
% Build normalization requant config from parameter map.
% Input:
%   task_param_map (containers.Map)
%       parameter key-value map
% Output:
%   cfg (struct)
%       normalization requant config
if ~isa(task_param_map, 'containers.Map')
    error('CNormalizationRequantServiceConfig:InvalidParamMap', ...
        'Expected task_param_map to be a containers.Map.');
end

base_cfg_struct = struct();
base_cfg_struct.spec_dir_path = CParamMapUtils.getRequired(task_param_map, CPTMdecoderWorkflowParamKeys.PARAM_SPEC_DIR_PATH, ...
    'spectrum directory path', 'CNormalizationRequantServiceConfig');

ms1_tol_value = CParamMapUtils.getRequiredNumber(task_param_map, CPTMdecoderWorkflowParamKeys.PARAM_MS1_TOLERANCE_VALUE, ...
    'MS1 tolerance value', 'CNormalizationRequantServiceConfig');
ms1_tol_type = strtrim(CParamMapUtils.getRequired(task_param_map, CPTMdecoderWorkflowParamKeys.PARAM_MS1_TOLERANCE_TYPE, ...
    'MS1 tolerance type', 'CNormalizationRequantServiceConfig'));
base_cfg_struct.ms1_tolerance = struct('value', ms1_tol_value, 'isppm', strcmpi(ms1_tol_type, 'PPM'));

base_cfg_struct.alpha = CParamMapUtils.getRequiredNumber(task_param_map, CPTMdecoderWorkflowParamKeys.PARAM_ALPHA, ...
    'noise filtering threshold', 'CNormalizationRequantServiceConfig');
base_cfg_struct.result_filter_threshold = CParamMapUtils.getRequiredNumber(task_param_map, CPTMdecoderWorkflowParamKeys.PARAM_RESULT_FILTER_THRESHOLD, ...
    'result filter threshold', 'CNormalizationRequantServiceConfig');
base_cfg_struct.output_dir_path = CParamMapUtils.getRequired(task_param_map, CPTMdecoderWorkflowParamKeys.PARAM_OUTPUT_DIR_PATH, ...
    'output directory', 'CNormalizationRequantServiceConfig');
base_cfg_struct.min_MSMS_num = CParamMapUtils.getOptionalNumber(task_param_map, CPTMdecoderWorkflowParamKeys.PARAM_MIN_MSMS_NUM, 1, 'CNormalizationRequantServiceConfig');
base_cfg_struct.checked_peptides_res_path = CParamMapUtils.getOptional(task_param_map, CPTMdecoderWorkflowParamKeys.PARAM_CHECKED_PEPTIDES_RES_PATH, []);

cfg = struct();
cfg.msms_cfg = base_cfg_struct;
cfg.output_file_name = CParamMapUtils.getOptional(task_param_map, ...
    CPTMdecoderWorkflowParamKeys.PARAM_NORM_REQUANT_OUTPUT_FILE_NAME, 'peptide4normalization_requant.txt');
end
