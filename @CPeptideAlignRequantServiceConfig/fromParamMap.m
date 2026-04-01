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

cfg.alpha = CParamMapUtils.getRequiredNumber(task_param_map, CPTMdecoderWorkflowParamKeys.PARAM_ALPHA, 'noise filtering threshold', 'CPeptideAlignRequantServiceConfig');
cfg.result_filter_threshold = CParamMapUtils.getRequiredNumber(task_param_map, CPTMdecoderWorkflowParamKeys.PARAM_RESULT_FILTER_THRESHOLD, 'result filter threshold', 'CPeptideAlignRequantServiceConfig');
cfg.fasta_file_path = CParamMapUtils.getRequired(task_param_map, CPTMdecoderWorkflowParamKeys.PARAM_FASTA_FILE_PATH, 'FASTA file path', 'CPeptideAlignRequantServiceConfig');
cfg.regular_express = CParamMapUtils.getRequired(task_param_map, CPTMdecoderWorkflowParamKeys.PARAM_REGULAR_EXPRESS, 'regular expression', 'CPeptideAlignRequantServiceConfig');
cfg.output_dir_path = CParamMapUtils.getOptional(task_param_map, CPTMdecoderWorkflowParamKeys.PARAM_OUTPUT_DIR_PATH, '', 'CPeptideAlignRequantServiceConfig');

cfg.filtered_res_file_path = CParamMapUtils.getOptional(task_param_map, CPTMdecoderWorkflowParamKeys.PARAM_FILTERED_RES_FILE_PATH, '', 'CPeptideAlignRequantServiceConfig');
cfg.msms_res_path = CParamMapUtils.getOptional(task_param_map, CPTMdecoderWorkflowParamKeys.PARAM_MSMS_RES_PATH, '', 'CPeptideAlignRequantServiceConfig');
cfg.peptide_quant_res_path = CParamMapUtils.getOptional(task_param_map, CPTMdecoderWorkflowParamKeys.PARAM_PEPTIDE_QUANT_RES_PATH, '', 'CPeptideAlignRequantServiceConfig');
cfg.min_MSMS_num = CParamMapUtils.getOptionalNumber(task_param_map, CPTMdecoderWorkflowParamKeys.PARAM_MIN_MSMS_NUM, 1, 'CPeptideAlignRequantServiceConfig');
cfg.alignment_report_path = CParamMapUtils.getOptional(task_param_map, ...
    CPTMdecoderWorkflowParamKeys.PARAM_ALIGN_REPORT_PATH, '', 'CPeptideAlignRequantServiceConfig');
cfg.requant_output_path = CParamMapUtils.getOptional(task_param_map, ...
    CPTMdecoderWorkflowParamKeys.PARAM_ALIGN_REQUANT_OUTPUT_PATH, '', 'CPeptideAlignRequantServiceConfig');
cfg.align_requant_rt_stats_path = CParamMapUtils.getOptional(task_param_map, ...
    CPTMdecoderWorkflowParamKeys.PARAM_ALIGN_REQUANT_RT_STATS_PATH, '', 'CPeptideAlignRequantServiceConfig');

cfg.align_strategy_obj = parseAlignStrategyFromMap(task_param_map);
cfg.align_options = parseAlignOptionsFromMap(task_param_map);
cfg.msms_stability_filter = CMsmsStabilityFilterConfig.fromParamMap(task_param_map, 'CPeptideAlignRequantServiceConfig');
cfg = CPeptideAlignRequantServiceConfig.finalize(cfg);
end


function align_strategy_obj = parseAlignStrategyFromMap(task_param_map)
% Build run-alignment strategy object from parameter map.
strategy_type = lower(strtrim(CParamMapUtils.getOptional(task_param_map, ...
    CPTMdecoderWorkflowParamKeys.PARAM_ALIGN_STRATEGY_TYPE, 'reference', 'CPeptideAlignRequantServiceConfig')));

switch strategy_type
    case 'reference'
        reference_raw = CParamMapUtils.getOptional(task_param_map, ...
            CPTMdecoderWorkflowParamKeys.PARAM_ALIGN_REFERENCE_RAW, '', 'CPeptideAlignRequantServiceConfig');
        align_strategy_obj = ReferenceRunAlignStrategy(reference_raw);
    case 'pairwise'
        pair_num = CParamMapUtils.getRequiredNumber(task_param_map, ...
            CPTMdecoderWorkflowParamKeys.PARAM_ALIGN_PAIR_NUM, ...
            'number of alignment pairs', 'CPeptideAlignRequantServiceConfig');
        if pair_num < 1 || floor(pair_num) ~= pair_num
            error('CPeptideAlignRequantServiceConfig:InvalidAlignPairNum', ...
                'align_pair_num must be a positive integer.');
        end
        pair_list = cell(pair_num, 2);
        for idx = 1:pair_num
            pair_key = [CPTMdecoderWorkflowParamKeys.PARAM_PREFIX_ALIGN_PAIR, num2str(idx)];
            pair_str = CParamMapUtils.getRequired(task_param_map, pair_key, ...
                'alignment pair in format left_raw|right_raw', 'CPeptideAlignRequantServiceConfig');
            split_pair = strsplit(pair_str, '|');
            if numel(split_pair) ~= 2
                error('CPeptideAlignRequantServiceConfig:InvalidAlignPair', ...
                    'Param ''%s'' must be in format left_raw|right_raw.', pair_key);
            end
            pair_list{idx, 1} = strtrim(split_pair{1});
            pair_list{idx, 2} = strtrim(split_pair{2});
        end
        align_strategy_obj = PairwiseRunAlignStrategy(pair_list);
    otherwise
        error('CPeptideAlignRequantServiceConfig:InvalidAlignStrategyType', ...
            'Unsupported align_strategy_type: %s (expected reference or pairwise).', strategy_type);
end
end


function align_options = parseAlignOptionsFromMap(task_param_map)
% Build align options struct from parameter map.
align_options = struct();

align_options.min_psm = CParamMapUtils.getOptionalNumber(task_param_map, ...
    CPTMdecoderWorkflowParamKeys.PARAM_ALIGN_MIN_PSM, [], 'CPeptideAlignRequantServiceConfig');
align_options.num_bins = CParamMapUtils.getOptionalNumber(task_param_map, ...
    CPTMdecoderWorkflowParamKeys.PARAM_ALIGN_NUM_BINS, [], 'CPeptideAlignRequantServiceConfig');
align_options.min_per_bin = CParamMapUtils.getOptionalNumber(task_param_map, ...
    CPTMdecoderWorkflowParamKeys.PARAM_ALIGN_MIN_PER_BIN, [], 'CPeptideAlignRequantServiceConfig');
align_options.outlier_k = CParamMapUtils.getOptionalNumber(task_param_map, ...
    CPTMdecoderWorkflowParamKeys.PARAM_ALIGN_OUTLIER_K, [], 'CPeptideAlignRequantServiceConfig');
align_options.outlier_method = CParamMapUtils.getOptional(task_param_map, ...
    CPTMdecoderWorkflowParamKeys.PARAM_ALIGN_OUTLIER_METHOD, [], 'CPeptideAlignRequantServiceConfig');
align_options.rt_sigma = CParamMapUtils.getOptionalNumber(task_param_map, ...
    CPTMdecoderWorkflowParamKeys.PARAM_ALIGN_RT_SIGMA, [], 'CPeptideAlignRequantServiceConfig');
align_options.max_rt_residual = CParamMapUtils.getOptionalNumber(task_param_map, ...
    CPTMdecoderWorkflowParamKeys.PARAM_ALIGN_MAX_RT_RESIDUAL, [], 'CPeptideAlignRequantServiceConfig');
align_options.dead_time_min = CParamMapUtils.getOptionalNumber(task_param_map, ...
    CPTMdecoderWorkflowParamKeys.PARAM_ALIGN_DEAD_TIME_MIN, [], 'CPeptideAlignRequantServiceConfig');
end
