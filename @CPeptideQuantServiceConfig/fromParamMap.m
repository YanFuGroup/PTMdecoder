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
cfg.msms_res_path = CParamMapUtils.getOptional(task_param_map, CPTMdecoderWorkflowParamKeys.PARAM_MSMS_RES_PATH, [], 'CPeptideQuantServiceConfig');
cfg.msms_stability_filter = parseMsmsStabilityFilter(task_param_map);
end


function filter_cfg = parseMsmsStabilityFilter(task_param_map)
% Build MSMS stability filter config for peptide-level quant stage.
% Input:
%   task_param_map (containers.Map)
%       parameter key-value map for peptide-level quant service
% Output:
%   filter_cfg (struct)
%       normalized filter config:
%       - enabled (logical)
%       - min_jaccard (1 x 1 double or [])
%       - min_support_frequency (1 x 1 double or [])
%       - max_abundance_mad (1 x 1 double or [])
%       - nan_as_fail (logical)
filter_cfg = struct();
filter_cfg.enabled = CParamMapUtils.getOptionalLogical(task_param_map, ...
    CPTMdecoderWorkflowParamKeys.PARAM_PEPTIDE_QUANT_STABILITY_FILTER_ON, ...
    false, 'CPeptideQuantServiceConfig');
filter_cfg.min_jaccard = CParamMapUtils.getOptionalNumber(task_param_map, ...
    CPTMdecoderWorkflowParamKeys.PARAM_PEPTIDE_QUANT_MIN_JACCARD_STABILITY, ...
    [], 'CPeptideQuantServiceConfig');
filter_cfg.min_support_frequency = CParamMapUtils.getOptionalNumber(task_param_map, ...
    CPTMdecoderWorkflowParamKeys.PARAM_PEPTIDE_QUANT_MIN_SUPPORT_FREQUENCY, ...
    [], 'CPeptideQuantServiceConfig');
filter_cfg.max_abundance_mad = CParamMapUtils.getOptionalNumber(task_param_map, ...
    CPTMdecoderWorkflowParamKeys.PARAM_PEPTIDE_QUANT_MAX_ABUNDANCE_MAD, ...
    [], 'CPeptideQuantServiceConfig');
filter_cfg.nan_as_fail = CParamMapUtils.getOptionalLogical(task_param_map, ...
    CPTMdecoderWorkflowParamKeys.PARAM_PEPTIDE_QUANT_STABILITY_NAN_AS_FAIL, ...
    true, 'CPeptideQuantServiceConfig');

if ~isempty(filter_cfg.min_jaccard) && (filter_cfg.min_jaccard < 0 || filter_cfg.min_jaccard > 1)
    CLogger.error('CPeptideQuantServiceConfig:InvalidMinJaccardStability', ...
        'peptide_quant_min_jaccard_stability must be within [0, 1].');
end
if ~isempty(filter_cfg.min_support_frequency) && ...
        (filter_cfg.min_support_frequency < 0 || filter_cfg.min_support_frequency > 1)
    CLogger.error('CPeptideQuantServiceConfig:InvalidMinSupportFrequency', ...
        'peptide_quant_min_support_frequency must be within [0, 1].');
end
if ~isempty(filter_cfg.max_abundance_mad) && filter_cfg.max_abundance_mad < 0
    CLogger.error('CPeptideQuantServiceConfig:InvalidMaxAbundanceMad', ...
        'peptide_quant_max_abundance_mad must be greater than or equal to 0.');
end
end
