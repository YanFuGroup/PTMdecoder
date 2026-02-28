function cfgs = fromParamMap(task_param_map)
% Build merge-to-pair configs from parameter map.
% Input:
%   task_param_map (containers.Map)
%       parameter key-value map for pair-level merge
% Output:
%   cfgs (cell)
%       cell array of CMergeEachPairConfig objects
if ~isa(task_param_map, 'containers.Map')
    error('CMergeEachPairConfig:InvalidParamMap', ...
        'Expected task_param_map to be a containers.Map.');
end

left_right_out_num = CParamMapUtils.getRequiredNumber(task_param_map, CPTMdecoderWorkflowParamKeys.PARAM_LEFT_RIGHT_OUT_NUM, 'number of pairwise comparisons');
left_name = CParamMapUtils.getRequired(task_param_map, CPTMdecoderWorkflowParamKeys.PARAM_LEFT_NAME, 'left group name');
right_name = CParamMapUtils.getRequired(task_param_map, CPTMdecoderWorkflowParamKeys.PARAM_RIGHT_NAME, 'right group name');
ignore_str = CParamMapUtils.getRequired(task_param_map, CPTMdecoderWorkflowParamKeys.PARAM_IGNORE_STRINGS_PAIR_LEVEL, 'ignore strings for pair-level merge');
ignore_strings = CParamMapUtils.parseQuotedList(ignore_str, ';');

cfgs = cell(left_right_out_num, 1);
for idx = 1:left_right_out_num
    key_name = [CPTMdecoderWorkflowParamKeys.PARAM_PREFIX_LEFT_RIGHT_OUT, num2str(idx)];
    left_right_out_str = CParamMapUtils.getRequired(task_param_map, key_name, 'pairwise input/output mapping');
    split_str = strsplit(left_right_out_str, {'|', '>'});
    pair_row = strtrim(split_str);
    cfgs{idx} = CMergeEachPairConfig.fromPairRow(pair_row, left_name, right_name, ignore_strings);
end
end
