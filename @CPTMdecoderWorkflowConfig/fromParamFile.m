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

cfg.stages = {};

% msms_action is an input intent flag, not a stage runtime action.
% Runtime dispatch only depends on explicit stage names.
msms_action = CPTMdecoderWorkflowConfig.resolveMsmsWorkflowActionFromMap(task_param_map);
if ~isempty(msms_action)
    if strcmp(msms_action, CPTMdecoderWorkflowConfig.ACTION_MSMS_PEPTIDE)
        msms_cfg = CMSMSLevelServiceConfig.fromParamMap(task_param_map);
        pep_quant_cfg = CPeptideQuantServiceConfig.fromParamMap(task_param_map);
        cfg.stages{end + 1} = CPTMdecoderWorkflowConfig.makeStage( ...
            CPTMdecoderWorkflowConfig.STAGE_MSMS_LEVEL, CPTMdecoderWorkflowConfig.ACTION_MSMS_ONLY, msms_cfg, true);
        cfg.stages{end + 1} = CPTMdecoderWorkflowConfig.makeStage( ...
            CPTMdecoderWorkflowConfig.STAGE_PEPTIDE_QUANT, CPTMdecoderWorkflowConfig.ACTION_PEPTIDE_ONLY, pep_quant_cfg, true);
    elseif strcmp(msms_action, CPTMdecoderWorkflowConfig.ACTION_MSMS_ONLY)
        msms_cfg = CMSMSLevelServiceConfig.fromParamMap(task_param_map);
        cfg.stages{end + 1} = CPTMdecoderWorkflowConfig.makeStage( ...
            CPTMdecoderWorkflowConfig.STAGE_MSMS_LEVEL, msms_action, msms_cfg, true);
    elseif strcmp(msms_action, CPTMdecoderWorkflowConfig.ACTION_PEPTIDE_ONLY)
        pep_quant_cfg = CPeptideQuantServiceConfig.fromParamMap(task_param_map);
        cfg.stages{end + 1} = CPTMdecoderWorkflowConfig.makeStage( ...
            CPTMdecoderWorkflowConfig.STAGE_PEPTIDE_QUANT, msms_action, pep_quant_cfg, true);
    elseif strcmp(msms_action, CPTMdecoderWorkflowConfig.ACTION_PEPTIDE_REQUANT)
        pep_requant_cfg = CPeptideRequantServiceConfig.fromParamMap(task_param_map);
        cfg.stages{end + 1} = CPTMdecoderWorkflowConfig.makeStage( ...
            CPTMdecoderWorkflowConfig.STAGE_PEPTIDE_REQUANT, msms_action, pep_requant_cfg, true);
    elseif strcmp(msms_action, CPTMdecoderWorkflowConfig.ACTION_PEPTIDE_ALIGN_REQUANT)
        pep_align_requant_cfg = CPeptideAlignRequantServiceConfig.fromParamMap(task_param_map);
        cfg.stages{end + 1} = CPTMdecoderWorkflowConfig.makeStage( ...
            CPTMdecoderWorkflowConfig.STAGE_PEPTIDE_ALIGN_REQUANT, msms_action, pep_align_requant_cfg, true);
    else
        error('CPTMdecoderWorkflowConfig:UnsupportedMsmsAction', ...
            'Unsupported msms action: %s', msms_action);
    end
end

if CPTMdecoderWorkflowConfig.getFlag(task_param_map, CPTMdecoderWorkflowParamKeys.PARAM_NORM_PEPTIDE_QUANT_ON)
    norm_quant_cfg = CNormalizationQuantServiceConfig.fromParamMap(task_param_map);
    cfg.stages{end + 1} = CPTMdecoderWorkflowConfig.makeStage( ...
        CPTMdecoderWorkflowConfig.STAGE_NORM_PEPTIDE_QUANT, CPTMdecoderWorkflowConfig.ACTION_NORM_PEPTIDE_QUANT, ...
        norm_quant_cfg, true);
end

if CPTMdecoderWorkflowConfig.getFlag(task_param_map, CPTMdecoderWorkflowParamKeys.PARAM_NORM_PEPTIDE_REQUANT_ON)
    norm_requant_cfg = CNormalizationRequantServiceConfig.fromParamMap(task_param_map);
    cfg.stages{end + 1} = CPTMdecoderWorkflowConfig.makeStage( ...
        CPTMdecoderWorkflowConfig.STAGE_NORM_PEPTIDE_REQUANT, CPTMdecoderWorkflowConfig.ACTION_NORM_PEPTIDE_REQUANT, ...
        norm_requant_cfg, true);
end

if CPTMdecoderWorkflowConfig.getFlag(task_param_map, CPTMdecoderWorkflowParamKeys.PARAM_SITE_LEVEL_ON)
    site_cfg = CSiteLevelPipelineConfig.fromParamMap(task_param_map);
    cfg.stages{end + 1} = CPTMdecoderWorkflowConfig.makeStage( ...
        CPTMdecoderWorkflowConfig.STAGE_SITE_LEVEL, CPTMdecoderWorkflowConfig.ACTION_SITE_SUMMARY, ...
        site_cfg, true);
end

if CPTMdecoderWorkflowConfig.getFlag(task_param_map, CPTMdecoderWorkflowParamKeys.PARAM_MERGE_TO_PAIR_LEVEL_ON)
    merge_to_pair_cfgs = CMergeEachPairConfig.fromParamMap(task_param_map);
    cfg.stages{end + 1} = CPTMdecoderWorkflowConfig.makeStage( ...
        CPTMdecoderWorkflowConfig.STAGE_MERGE_TO_PAIR_LEVEL, CPTMdecoderWorkflowConfig.ACTION_MERGE_EACH_PAIR, ...
        merge_to_pair_cfgs, true);
end

if CPTMdecoderWorkflowConfig.getFlag(task_param_map, CPTMdecoderWorkflowParamKeys.PARAM_MERGE_PAIRS_LEVEL_ON)
    merge_pairs_cfg = CMergePairsConfig.fromParamMap(task_param_map);
    cfg.stages{end + 1} = CPTMdecoderWorkflowConfig.makeStage( ...
        CPTMdecoderWorkflowConfig.STAGE_MERGE_PAIRS_LEVEL, CPTMdecoderWorkflowConfig.ACTION_MERGE_PAIRS, ...
        merge_pairs_cfg, true);
end

obj = CPTMdecoderWorkflowConfig(cfg);
end