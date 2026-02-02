function drawGather(obj, pep_rtrange_map, dir_save, color_map, legend_map)
% Draw the XIC for gathered peptides using manually-checked rt range
% Input:
%   obj (CIMPGatherQuant)
%       Quantification aggregator instance
%   pep_rtrange_map (containers.Map)
%       map of [modified peptide _ charge _ raw file name] -> [rt_start, rt_end, check_label]
%   dir_save (1 x 1 char/string)
%       directory to save the figures
%   color_map (containers.Map or [])
%       color map (key: imp name, value: RGB 1x3)
%   legend_map (containers.Map or [])
%       legend map (key: imp name, value: display string)

% Check the input arguments
if nargin < 5
    legend_map = [];
end
if nargin < 4
    color_map = [];
end


state = struct('dir_save', dir_save, 'color_map', color_map, 'legend_map', legend_map);

% Do the same operation for gathered PSM for every raw file
iterate_imp_groups(obj, pep_rtrange_map, @handle_group_charge, state);
end



function state = handle_group_charge(state, obj, raw_name, group_imp_name, group_ratio, group_rts, group_inten, ...
        ~, low_mz_bound, high_mz_bound, selected_charge, charge_group_idxs, current_imp_rt_range)
if isempty(current_imp_rt_range) || all(cellfun(@isempty,current_imp_rt_range))
    % All of the IMPs are removed in manual checking
    return;
end

% Draw for this group
obj.draw_each_group(raw_name,...
    group_ratio(charge_group_idxs,:),...
    group_rts(charge_group_idxs,:),...
    group_inten(charge_group_idxs,:),...
    low_mz_bound,high_mz_bound, ...
    selected_charge,current_imp_rt_range,...
    group_imp_name, state.dir_save, state.color_map, state.legend_map);
end