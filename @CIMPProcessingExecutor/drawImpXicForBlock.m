function drawImpXicForBlock(obj, rawIdentManager, pep_rtrange_map, dir_save, color_map, legend_map, xic_layout)
% Draw XICs for IMP groups using checked RT ranges
% Input:
%   obj (CIMPProcessingExecutor)
%       processing executor instance
%   rawIdentManager (CIMPRawIdentManager)
%       per-raw identification store manager
%   pep_rtrange_map (containers.Map)
%       map of [modified peptide _ charge _ raw file name] -> [rt_start, rt_end, check_label]
%   dir_save (1 x 1 char/string)
%       directory to save the figures
%   color_map (containers.Map or [])
%       color map (key: imp name, value: RGB 1x3)
%   legend_map (containers.Map or [])
%       legend map (key: imp name, value: display string)
%   xic_layout (struct or [], optional)
%       optional XIC layout override

if nargin < 7
    xic_layout = [];
end

state = struct('dir_save', dir_save, 'color_map', color_map, ...
    'legend_map', legend_map, 'xic_layout', xic_layout);
[raw_names, raw_ident_stores] = rawIdentManager.getEntries();

groupAggregator = CIMPGroupAggregator(obj.m_ms1_tolerance);
groupAggregator.aggregate(raw_names, raw_ident_stores, pep_rtrange_map, ...
    @(state_in, group) obj.onGroupDrawXic(state_in, group), state);
end
