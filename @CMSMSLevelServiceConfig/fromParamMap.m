function cfg = fromParamMap(task_param_map)
% Build MSMS-level service config struct from parameter map.
% Input:
%   task_param_map (containers.Map)
%       parameter key-value map for MSMS-level service
% Output:
%   cfg (struct)
%       normalized service config
if ~isa(task_param_map, 'containers.Map')
    error('CMSMSLevelServiceConfig:InvalidParamMap', ...
        'Expected task_param_map to be a containers.Map.');
end

cfg = struct();
cfg.mod_file_path = CParamMapUtils.getRequired(task_param_map, CPTMdecoderWorkflowParamKeys.PARAM_MOD_FILE_PATH, 'modification file path', 'CMSMSLevelServiceConfig');
cfg.fixed_mod = CParamMapUtils.getRequired(task_param_map, CPTMdecoderWorkflowParamKeys.PARAM_FIXED_MOD, 'fixed modifications', 'CMSMSLevelServiceConfig');
cfg.variable_mod = CParamMapUtils.getRequired(task_param_map, CPTMdecoderWorkflowParamKeys.PARAM_VARIABLE_MOD, 'variable modifications', 'CMSMSLevelServiceConfig');
cfg.spec_dir_path = CParamMapUtils.getRequired(task_param_map, CPTMdecoderWorkflowParamKeys.PARAM_SPEC_DIR_PATH, 'spectrum directory path', 'CMSMSLevelServiceConfig');

ms1_tol_value = CParamMapUtils.getRequiredNumber(task_param_map, CPTMdecoderWorkflowParamKeys.PARAM_MS1_TOLERANCE_VALUE, 'MS1 tolerance value', 'CMSMSLevelServiceConfig');
ms1_tol_type = strtrim(CParamMapUtils.getRequired(task_param_map, CPTMdecoderWorkflowParamKeys.PARAM_MS1_TOLERANCE_TYPE, 'MS1 tolerance type', 'CMSMSLevelServiceConfig'));
cfg.ms1_tolerance = struct('value', ms1_tol_value, 'isppm', strcmpi(ms1_tol_type, 'PPM'));

cfg.ms2_tolerance = CParamMapUtils.getRequiredNumber(task_param_map, CPTMdecoderWorkflowParamKeys.PARAM_MS2_TOLERANCE, 'MS2 tolerance value', 'CMSMSLevelServiceConfig');
cfg.alpha = CParamMapUtils.getRequiredNumber(task_param_map, CPTMdecoderWorkflowParamKeys.PARAM_ALPHA, 'noise filtering threshold', 'CMSMSLevelServiceConfig');
cfg.fasta_file_path = CParamMapUtils.getRequired(task_param_map, CPTMdecoderWorkflowParamKeys.PARAM_FASTA_FILE_PATH, 'FASTA file path', 'CMSMSLevelServiceConfig');
cfg.regular_express = CParamMapUtils.getRequired(task_param_map, CPTMdecoderWorkflowParamKeys.PARAM_REGULAR_EXPRESS, 'regular expression', 'CMSMSLevelServiceConfig');
cfg.filtered_res_file_path = CParamMapUtils.getOptional(task_param_map, CPTMdecoderWorkflowParamKeys.PARAM_FILTERED_RES_FILE_PATH, '', 'CMSMSLevelServiceConfig');
cfg.model = CParamMapUtils.getRequiredNumber(task_param_map, CPTMdecoderWorkflowParamKeys.PARAM_MODEL, 'quantification model', 'CMSMSLevelServiceConfig');
cfg.method = CParamMapUtils.getRequiredNumber(task_param_map, CPTMdecoderWorkflowParamKeys.PARAM_METHOD, 'solver method', 'CMSMSLevelServiceConfig');

if task_param_map.isKey(CPTMdecoderWorkflowParamKeys.PARAM_LAMBDA)
    cfg.lambda = CParamMapUtils.getRequiredNumber(task_param_map, CPTMdecoderWorkflowParamKeys.PARAM_LAMBDA, 'lasso parameter lambda', 'CMSMSLevelServiceConfig');
else
    cfg.lambda = 0;
end
if cfg.method == 2 && ~task_param_map.isKey(CPTMdecoderWorkflowParamKeys.PARAM_LAMBDA)
    error('CMSMSLevelServiceConfig:MissingLambda', ...
        'The lasso parameter ''lambda'' is required when method=2.');
end

cfg.result_filter_threshold = CParamMapUtils.getRequiredNumber(task_param_map, CPTMdecoderWorkflowParamKeys.PARAM_RESULT_FILTER_THRESHOLD, 'result filter threshold', 'CMSMSLevelServiceConfig');
cfg.enzyme_name = CParamMapUtils.getRequired(task_param_map, CPTMdecoderWorkflowParamKeys.PARAM_ENZYME_NAME, 'enzyme name', 'CMSMSLevelServiceConfig');
cfg.enzyme_limits = str2num(CParamMapUtils.getRequired(task_param_map, CPTMdecoderWorkflowParamKeys.PARAM_ENZYME_LIMIT_C_TERM_POSSIBLE_MOD, 'enzyme limit C-term possible modifications', 'CMSMSLevelServiceConfig')); %#ok<ST2NM>
cfg.output_dir_path = CParamMapUtils.getRequired(task_param_map, CPTMdecoderWorkflowParamKeys.PARAM_OUTPUT_DIR_PATH, 'output directory path', 'CMSMSLevelServiceConfig');
cfg.msms_report_vif_on = CParamMapUtils.getOptionalNumber(task_param_map, CPTMdecoderWorkflowParamKeys.PARAM_MSMS_REPORT_VIF_ON, 0, 'CMSMSLevelServiceConfig');
if ~(cfg.msms_report_vif_on == 0 || cfg.msms_report_vif_on == 1)
    error('CMSMSLevelServiceConfig:InvalidMsmsReportVifOn', ...
        'Param ''msms_report_vif_on'' must be 0 or 1.');
end

if task_param_map.isKey(CPTMdecoderWorkflowParamKeys.PARAM_N_RESAMPLES)
    num_resamples = CParamMapUtils.getRequiredNumber(task_param_map, CPTMdecoderWorkflowParamKeys.PARAM_N_RESAMPLES, 'number of stability resamples', 'CMSMSLevelServiceConfig');
    if num_resamples <= 0 || num_resamples ~= floor(num_resamples)
        error('CMSMSLevelServiceConfig:InvalidNResamples', ...
            'Param ''n_resamples'' must be a positive integer.');
    end

    random_seed = CParamMapUtils.getOptionalNumber(task_param_map, CPTMdecoderWorkflowParamKeys.PARAM_RANDOM_SEED, 42, 'CMSMSLevelServiceConfig');
    use_parallel = CParamMapUtils.getOptionalLogical( ...
        task_param_map, ...
        CPTMdecoderWorkflowParamKeys.PARAM_STABILITY_PARALLEL_ON, ...
        false, ...
        'CMSMSLevelServiceConfig');

    cfg.stability_options = struct( ...
        'n_resamples', num_resamples, ...
        'random_seed', random_seed, ...
        'use_parallel', use_parallel, ...
        'relative_threshold', cfg.result_filter_threshold);
end

cfg.min_MSMS_num = CParamMapUtils.getOptionalNumber(task_param_map, CPTMdecoderWorkflowParamKeys.PARAM_MIN_MSMS_NUM, 1, 'CMSMSLevelServiceConfig');
cfg.min_peptide_length = CParamMapUtils.getOptionalNumber(task_param_map, CPTMdecoderWorkflowParamKeys.PARAM_PEPTIDE_MIN_LENGTH, 7, 'CMSMSLevelServiceConfig');
cfg.max_peptide_length = CParamMapUtils.getOptionalNumber(task_param_map, CPTMdecoderWorkflowParamKeys.PARAM_PEPTIDE_MAX_LENGTH, 40, 'CMSMSLevelServiceConfig');
cfg.max_mod_per_peptide = CParamMapUtils.getOptionalNumber(task_param_map, CPTMdecoderWorkflowParamKeys.PARAM_MAX_MOD_PER_PEPTIDE, 5, 'CMSMSLevelServiceConfig');

if cfg.min_peptide_length <= 0 || cfg.min_peptide_length ~= floor(cfg.min_peptide_length)
    error('CMSMSLevelServiceConfig:InvalidMinPeptideLength', ...
        'Param ''peptide_min_length'' must be a positive integer.');
end
if cfg.max_peptide_length <= 0 || cfg.max_peptide_length ~= floor(cfg.max_peptide_length)
    error('CMSMSLevelServiceConfig:InvalidMaxPeptideLength', ...
        'Param ''peptide_max_length'' must be a positive integer.');
end
if cfg.min_peptide_length > cfg.max_peptide_length
    error('CMSMSLevelServiceConfig:InvalidPeptideLengthRange', ...
        'Param ''peptide_min_length'' must be less than or equal to ''peptide_max_length''.');
end
if cfg.max_mod_per_peptide <= 0 || cfg.max_mod_per_peptide ~= floor(cfg.max_mod_per_peptide)
    error('CMSMSLevelServiceConfig:InvalidMaxModPerPeptide', ...
        'Param ''max_mod_per_peptide'' must be a positive integer.');
end

cfg.ion_types = [1,2];
cfg.case_penalty_intens = 'intens_sum';
cfg.grid_penalty_intens = 'intens_sum';
cfg.case_OLS_intens_weight = 'none';
cfg.pep_spec_file_path = CParamMapUtils.getRequired(task_param_map, CPTMdecoderWorkflowParamKeys.PARAM_PEP_SPEC_FILE_PATH, 'peptide-spectrum file path', 'CMSMSLevelServiceConfig');
cfg.checked_peptides_res_path = [];
end
