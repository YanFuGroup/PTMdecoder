function cfg = fromParamMap(task_param_map)
% Build merge-pairs config from parameter map.
% Input:
%   task_param_map (containers.Map)
%       parameter key-value map for final merge
% Output:
%   cfg (CMergePairsConfig)
if ~isa(task_param_map, 'containers.Map')
    error('CMergePairsConfig:InvalidParamMap', ...
        'Expected task_param_map to be a containers.Map.');
end

pair_num = CParamMapUtils.getRequiredNumber(task_param_map, CPTMdecoderWorkflowParamKeys.PARAM_PAIR_NUM, 'number of pairs to merge');
result_paths = cell(pair_num, 1);
group_titles = cell(pair_num, 2);
for idx = 1:pair_num
    pair_key = [CPTMdecoderWorkflowParamKeys.PARAM_PREFIX_PAIR, num2str(idx)];
    lr_key = [CPTMdecoderWorkflowParamKeys.PARAM_PREFIX_LEFT_RIGHT_NAME, num2str(idx)];
    result_paths{idx} = CParamMapUtils.getRequired(task_param_map, pair_key, 'pair-level result path');

    lr_str = CParamMapUtils.getRequired(task_param_map, lr_key, 'left-right group names');
    split_lr = strsplit(lr_str, '|');
    group_titles{idx, 1} = strtrim(split_lr{1});
    group_titles{idx, 2} = strtrim(split_lr{2});
end

cfg_struct = struct();
cfg_struct.result_paths = result_paths;
cfg_struct.output_path = CParamMapUtils.getRequired(task_param_map, CPTMdecoderWorkflowParamKeys.PARAM_FINAL_OUTPUT_PATH, 'final merged output path');
cfg_struct.group_titles = group_titles;

cfg = CMergePairsConfig(cfg_struct);
end
