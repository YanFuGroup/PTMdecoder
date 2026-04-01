function obj = fromParamFile(param_file_path)
% Build workflow config directly from parameter file.
% Inputs:
%   param_file_path: (1 x 1 char)
%       Path to the parameter file containing all necessary parameters to build the workflow config.
% Outputs:
%   obj: (1 x 1 CPTMdecoderWorkflowConfig)
task_param_map = CWorkflowParamParser.parseFileToMap(param_file_path);

cfg = struct();
cfg.param_file_path = param_file_path;

logger_cfg = buildLoggerConfigFromParamMap(task_param_map);
CLogger.configure(logger_cfg);
CLogger.flush();

cfg.stages = {};

msms_flags = CPTMdecoderWorkflowConfig.resolveMsmsWorkflowFlagsFromMap(task_param_map);
if msms_flags.run_msms_level
    msms_cfg = CMSMSLevelServiceConfig.fromParamMap(task_param_map);
    cfg.stages{end + 1} = CPTMdecoderWorkflowConfig.makeStage( ...
        CPTMdecoderWorkflowConfig.STAGE_MSMS_LEVEL, msms_cfg, true);
end

if msms_flags.run_peptide_quant
    pep_quant_cfg = CPeptideQuantServiceConfig.fromParamMap(task_param_map);
    cfg.stages{end + 1} = CPTMdecoderWorkflowConfig.makeStage( ...
        CPTMdecoderWorkflowConfig.STAGE_PEPTIDE_QUANT, pep_quant_cfg, true);
end

if CPTMdecoderWorkflowConfig.getFlag(task_param_map, CPTMdecoderWorkflowParamKeys.PARAM_PEPTIDE_REQUANT_ON)
    pep_requant_cfg = CPeptideRequantServiceConfig.fromParamMap(task_param_map);
    cfg.stages{end + 1} = CPTMdecoderWorkflowConfig.makeStage( ...
        CPTMdecoderWorkflowConfig.STAGE_PEPTIDE_REQUANT, pep_requant_cfg, true);
end

if CPTMdecoderWorkflowConfig.getFlag(task_param_map, CPTMdecoderWorkflowParamKeys.PARAM_PEPTIDE_ALIGN_REQUANT_ON)
    pep_align_requant_cfg = CPeptideAlignRequantServiceConfig.fromParamMap(task_param_map);
    cfg.stages{end + 1} = CPTMdecoderWorkflowConfig.makeStage( ...
        CPTMdecoderWorkflowConfig.STAGE_PEPTIDE_ALIGN_REQUANT, pep_align_requant_cfg, true);
end

if CPTMdecoderWorkflowConfig.getFlag(task_param_map, CPTMdecoderWorkflowParamKeys.PARAM_NORM_PEPTIDE_QUANT_ON)
    norm_quant_cfg = CNormalizationQuantServiceConfig.fromParamMap(task_param_map);
    cfg.stages{end + 1} = CPTMdecoderWorkflowConfig.makeStage( ...
        CPTMdecoderWorkflowConfig.STAGE_NORM_PEPTIDE_QUANT, norm_quant_cfg, true);
end

if CPTMdecoderWorkflowConfig.getFlag(task_param_map, CPTMdecoderWorkflowParamKeys.PARAM_NORM_PEPTIDE_REQUANT_ON)
    norm_requant_cfg = CNormalizationRequantServiceConfig.fromParamMap(task_param_map);
    cfg.stages{end + 1} = CPTMdecoderWorkflowConfig.makeStage( ...
        CPTMdecoderWorkflowConfig.STAGE_NORM_PEPTIDE_REQUANT, norm_requant_cfg, true);
end

if CPTMdecoderWorkflowConfig.getFlag(task_param_map, CPTMdecoderWorkflowParamKeys.PARAM_SITE_LEVEL_ON)
    site_cfg = CSiteLevelPipelineConfig.fromParamMap(task_param_map);
    cfg.stages{end + 1} = CPTMdecoderWorkflowConfig.makeStage( ...
        CPTMdecoderWorkflowConfig.STAGE_SITE_LEVEL, site_cfg, true);
end

if CPTMdecoderWorkflowConfig.getFlag(task_param_map, CPTMdecoderWorkflowParamKeys.PARAM_SITE_LEVEL_DATASET_ON)
    site_dataset_cfg = CSiteLevelDatasetPipelineConfig.fromParamMap(task_param_map);
    cfg.stages{end + 1} = CPTMdecoderWorkflowConfig.makeStage( ...
        CPTMdecoderWorkflowConfig.STAGE_SITE_LEVEL_DATASET, site_dataset_cfg, true);
end

if CPTMdecoderWorkflowConfig.getFlag(task_param_map, CPTMdecoderWorkflowParamKeys.PARAM_MERGE_TO_PAIR_LEVEL_ON)
    merge_to_pair_cfgs = CMergeEachPairConfig.fromParamMap(task_param_map);
    cfg.stages{end + 1} = CPTMdecoderWorkflowConfig.makeStage( ...
        CPTMdecoderWorkflowConfig.STAGE_MERGE_TO_PAIR_LEVEL, merge_to_pair_cfgs, true);
end

if CPTMdecoderWorkflowConfig.getFlag(task_param_map, CPTMdecoderWorkflowParamKeys.PARAM_MERGE_PAIRS_LEVEL_ON)
    merge_pairs_cfg = CMergePairsConfig.fromParamMap(task_param_map);
    cfg.stages{end + 1} = CPTMdecoderWorkflowConfig.makeStage( ...
        CPTMdecoderWorkflowConfig.STAGE_MERGE_PAIRS_LEVEL, merge_pairs_cfg, true);
end

obj = CPTMdecoderWorkflowConfig(cfg);
end


function logger_cfg = buildLoggerConfigFromParamMap(task_param_map)
% Build logger configuration struct from the parameter map.
% Inputs:
%   task_param_map (containers.Map)
%       parameter key-value map
% Output:
%   logger_cfg (struct)
%       logger config struct with fields:
%       - enabled (logical)
%       - file_path (char)
%       - file_level (char)
%       - to_console (logical)
%       - console_level (char)

timestamp = datestr(now, 'yyyymmdd_HHMMSS');
log_file_name = sprintf('ptmdecoder_%s.log', timestamp);

logger_cfg = struct();
logger_cfg.enabled = CParamMapUtils.getOptionalLogical(task_param_map, ...
    CPTMdecoderWorkflowParamKeys.PARAM_LOG_ENABLED, true, 'CPTMdecoderWorkflowConfig');
configured_log_dir = CParamMapUtils.getOptional(task_param_map, ...
    CPTMdecoderWorkflowParamKeys.PARAM_LOG_FILE_DIR, '', 'CPTMdecoderWorkflowConfig');
if isempty(configured_log_dir)
    logger_cfg.file_path = fullfile(pwd, log_file_name);
else
    logger_cfg.file_path = fullfile(char(configured_log_dir), log_file_name);
end
logger_cfg.file_level = upper(strtrim(char(CParamMapUtils.getOptional(task_param_map, ...
    CPTMdecoderWorkflowParamKeys.PARAM_LOG_FILE_LEVEL, 'DEBUG', 'CPTMdecoderWorkflowConfig'))));
logger_cfg.to_console = CParamMapUtils.getOptionalLogical(task_param_map, ...
    CPTMdecoderWorkflowParamKeys.PARAM_LOG_TO_CONSOLE, true, 'CPTMdecoderWorkflowConfig');
logger_cfg.console_level = upper(strtrim(char(CParamMapUtils.getOptional(task_param_map, ...
    CPTMdecoderWorkflowParamKeys.PARAM_LOG_CONSOLE_LEVEL, 'INFO', 'CPTMdecoderWorkflowConfig'))));
end