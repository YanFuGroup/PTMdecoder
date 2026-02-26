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
    cfg_struct = CPTMdecoderWorkflowConfig.buildMsmsPepConfigStructFromMap(task_param_map, msms_action);
    if strcmp(msms_action, CPTMdecoderWorkflowConfig.ACTION_MSMS_PEPTIDE)
        msms_cfg = cfg_struct;
        cfg.stages{end + 1} = CPTMdecoderWorkflowConfig.makeStage( ...
            CPTMdecoderWorkflowConfig.STAGE_MSMS_LEVEL, CPTMdecoderWorkflowConfig.ACTION_MSMS_ONLY, msms_cfg, true);
        cfg.stages{end + 1} = CPTMdecoderWorkflowConfig.makeStage( ...
            CPTMdecoderWorkflowConfig.STAGE_PEPTIDE_QUANT, CPTMdecoderWorkflowConfig.ACTION_PEPTIDE_ONLY, msms_cfg, true);
    elseif strcmp(msms_action, CPTMdecoderWorkflowConfig.ACTION_MSMS_ONLY)
        cfg.stages{end + 1} = CPTMdecoderWorkflowConfig.makeStage( ...
            CPTMdecoderWorkflowConfig.STAGE_MSMS_LEVEL, msms_action, cfg_struct, true);
    elseif strcmp(msms_action, CPTMdecoderWorkflowConfig.ACTION_PEPTIDE_ONLY)
        cfg.stages{end + 1} = CPTMdecoderWorkflowConfig.makeStage( ...
            CPTMdecoderWorkflowConfig.STAGE_PEPTIDE_QUANT, msms_action, cfg_struct, true);
    elseif strcmp(msms_action, CPTMdecoderWorkflowConfig.ACTION_PEPTIDE_REQUANT)
        cfg.stages{end + 1} = CPTMdecoderWorkflowConfig.makeStage( ...
            CPTMdecoderWorkflowConfig.STAGE_PEPTIDE_REQUANT, msms_action, cfg_struct, true);
    else
        error('CPTMdecoderWorkflowConfig:UnsupportedMsmsAction', ...
            'Unsupported msms action: %s', msms_action);
    end
end

if CPTMdecoderWorkflowConfig.getFlag(task_param_map, CPTMdecoderWorkflowConfig.PARAM_NORM_PEPTIDE_QUANT_ON)
    norm_quant_cfg = CPTMdecoderWorkflowConfig.buildNormalizationQuantConfigFromMap(task_param_map);
    cfg.stages{end + 1} = CPTMdecoderWorkflowConfig.makeStage( ...
        CPTMdecoderWorkflowConfig.STAGE_NORM_PEPTIDE_QUANT, CPTMdecoderWorkflowConfig.ACTION_NORM_PEPTIDE_QUANT, ...
        norm_quant_cfg, true);
end

if CPTMdecoderWorkflowConfig.getFlag(task_param_map, CPTMdecoderWorkflowConfig.PARAM_NORM_PEPTIDE_REQUANT_ON)
    norm_requant_cfg = CPTMdecoderWorkflowConfig.buildNormalizationRequantConfigFromMap(task_param_map);
    cfg.stages{end + 1} = CPTMdecoderWorkflowConfig.makeStage( ...
        CPTMdecoderWorkflowConfig.STAGE_NORM_PEPTIDE_REQUANT, CPTMdecoderWorkflowConfig.ACTION_NORM_PEPTIDE_REQUANT, ...
        norm_requant_cfg, true);
end

if CPTMdecoderWorkflowConfig.getFlag(task_param_map, CPTMdecoderWorkflowConfig.PARAM_SITE_LEVEL_ON)
    site_cfg_struct = CPTMdecoderWorkflowConfig.buildSiteConfigStructFromMap(task_param_map);
    cfg.stages{end + 1} = CPTMdecoderWorkflowConfig.makeStage( ...
        CPTMdecoderWorkflowConfig.STAGE_SITE_LEVEL, CPTMdecoderWorkflowConfig.ACTION_SITE_SUMMARY, ...
        CSiteLevelPipelineConfig(site_cfg_struct), true);
end

if CPTMdecoderWorkflowConfig.getFlag(task_param_map, CPTMdecoderWorkflowConfig.PARAM_MERGE_TO_PAIR_LEVEL_ON)
    merge_to_pair_cfgs = CPTMdecoderWorkflowConfig.buildMergeToPairConfigsFromMap(task_param_map);
    cfg.stages{end + 1} = CPTMdecoderWorkflowConfig.makeStage( ...
        CPTMdecoderWorkflowConfig.STAGE_MERGE_TO_PAIR_LEVEL, CPTMdecoderWorkflowConfig.ACTION_MERGE_EACH_PAIR, ...
        merge_to_pair_cfgs, true);
end

if CPTMdecoderWorkflowConfig.getFlag(task_param_map, CPTMdecoderWorkflowConfig.PARAM_MERGE_PAIRS_LEVEL_ON)
    merge_pairs_cfg_struct = CPTMdecoderWorkflowConfig.buildMergePairsConfigStructFromMap(task_param_map);
    cfg.stages{end + 1} = CPTMdecoderWorkflowConfig.makeStage( ...
        CPTMdecoderWorkflowConfig.STAGE_MERGE_PAIRS_LEVEL, CPTMdecoderWorkflowConfig.ACTION_MERGE_PAIRS, ...
        CMergePairsConfig(merge_pairs_cfg_struct), true);
end

obj = CPTMdecoderWorkflowConfig(cfg);
end