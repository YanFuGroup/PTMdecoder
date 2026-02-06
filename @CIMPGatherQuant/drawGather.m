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
[raw_names, raw_ident_stores] = obj.getRawEntries();
for idx_raw = 1:numel(raw_names)
    state = obj.m_groupAggregator.aggregate(raw_names{idx_raw}, raw_ident_stores{idx_raw}, pep_rtrange_map, ...
        @(state_in, group) handle_group_charge(state_in, obj, group), state);
end
end



function state = handle_group_charge(state, obj, group)
if isempty(group.impRtRanges) || all(cellfun(@isempty, group.impRtRanges))
    % All of the IMPs are removed in manual checking
    return;
end

% Draw for this group
CIMPQuantPlotter.drawGroup(obj.m_cMs12DatasetIO, obj.m_minMSMSnum, group.rawName, ...
    group.ratio(group.chargeGroupIdxs,:), ...
    group.rts(group.chargeGroupIdxs,:), ...
    group.intensity(group.chargeGroupIdxs,:), ...
    group.lowMzBound, group.highMzBound, ...
    group.selectedCharge, group.impRtRanges, ...
    group.impNames, state.dir_save, state.color_map, state.legend_map);
end