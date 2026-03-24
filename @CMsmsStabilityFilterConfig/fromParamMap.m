function filter_cfg = fromParamMap(task_param_map, err_id_prefix)
% Build MSMS stability filter config from parameter map.
% Inputs:
%   task_param_map (containers.Map)
%       Parameter key-value map.
%   err_id_prefix (1 x 1 char/string)
%       Error id prefix used in CLogger.error tags.
% Outputs:
%   filter_cfg (1 x 1 struct)
%       Normalized filter config:
%       - enabled (logical)
%       - min_jaccard (1 x 1 double or [])
%       - min_support_frequency (1 x 1 double or [])
%       - max_abundance_mad (1 x 1 double or [])
%       - nan_as_fail (logical)
if nargin < 2 || isempty(err_id_prefix)
    err_id_prefix = 'CMsmsStabilityFilterConfig';
end

if ~isa(task_param_map, 'containers.Map')
    error([char(err_id_prefix), ':InvalidParamMap'], ...
        'Expected task_param_map to be a containers.Map.');
end

filter_cfg = struct();
filter_cfg.enabled = CParamMapUtils.getOptionalLogical(task_param_map, ...
    CPTMdecoderWorkflowParamKeys.PARAM_PEPTIDE_QUANT_STABILITY_FILTER_ON, ...
    false, err_id_prefix);
filter_cfg.min_jaccard = CParamMapUtils.getOptionalNumber(task_param_map, ...
    CPTMdecoderWorkflowParamKeys.PARAM_PEPTIDE_QUANT_MIN_JACCARD_STABILITY, ...
    [], err_id_prefix);
filter_cfg.min_support_frequency = CParamMapUtils.getOptionalNumber(task_param_map, ...
    CPTMdecoderWorkflowParamKeys.PARAM_PEPTIDE_QUANT_MIN_SUPPORT_FREQUENCY, ...
    [], err_id_prefix);
filter_cfg.max_abundance_mad = CParamMapUtils.getOptionalNumber(task_param_map, ...
    CPTMdecoderWorkflowParamKeys.PARAM_PEPTIDE_QUANT_MAX_ABUNDANCE_MAD, ...
    [], err_id_prefix);
filter_cfg.nan_as_fail = CParamMapUtils.getOptionalLogical(task_param_map, ...
    CPTMdecoderWorkflowParamKeys.PARAM_PEPTIDE_QUANT_STABILITY_NAN_AS_FAIL, ...
    true, err_id_prefix);

if ~isempty(filter_cfg.min_jaccard) && (filter_cfg.min_jaccard < 0 || filter_cfg.min_jaccard > 1)
    CLogger.error([char(err_id_prefix), ':InvalidMinJaccardStability'], ...
        'peptide_quant_min_jaccard_stability must be within [0, 1].');
end
if ~isempty(filter_cfg.min_support_frequency) && ...
        (filter_cfg.min_support_frequency < 0 || filter_cfg.min_support_frequency > 1)
    CLogger.error([char(err_id_prefix), ':InvalidMinSupportFrequency'], ...
        'peptide_quant_min_support_frequency must be within [0, 1].');
end
if ~isempty(filter_cfg.max_abundance_mad) && filter_cfg.max_abundance_mad < 0
    CLogger.error([char(err_id_prefix), ':InvalidMaxAbundanceMad'], ...
        'peptide_quant_max_abundance_mad must be greater than or equal to 0.');
end
end
