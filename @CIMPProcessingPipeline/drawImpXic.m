function drawImpXic(obj, rawManager, pep_rtrange_map, dir_save, color_map, legend_map)
% Draw XICs for IMP groups using checked RT ranges
% Input:
%   obj (CIMPProcessingPipeline)
%       processing pipeline instance
%   rawManager (CIMPRawIdentManager)
%       per-raw identification store manager
%   pep_rtrange_map (containers.Map)
%       map of [modified peptide _ charge _ raw file name] -> [rt_start, rt_end, check_label]
%   dir_save (1 x 1 char/string)
%       directory to save the figures
%   color_map (containers.Map or [])
%       color map (key: imp name, value: RGB 1x3)
%   legend_map (containers.Map or [])
%       legend map (key: imp name, value: display string)

state = struct('dir_save', dir_save, 'color_map', color_map, 'legend_map', legend_map);
[raw_names, raw_ident_stores] = rawManager.getEntries();

groupAggregator = CIMPGroupAggregator(obj.m_ms1_tolerance);
groupAggregator.aggregate(raw_names, raw_ident_stores, pep_rtrange_map, ...
    @(state_in, group) obj.drawGroupXic(state_in, group), state);
end
