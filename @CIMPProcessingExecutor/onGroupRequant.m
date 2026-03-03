function imp_records = onGroupRequant(obj, imp_records, group)
% Append re-quantification records for a single group
% Input:
%   obj (CIMPProcessingExecutor)
%       processing executor instance
%   imp_records (CIMPQuantRecord array)
%       accumulator for IMP records
%   group (CIMPGroup)
%       group payload for one raw and charge
% Output:
%   imp_records (CIMPQuantRecord array)
%       updated accumulator

if isempty(group.impRtRanges) || all(cellfun(@isempty, group.impRtRanges))
    return;
end

[has_nonzero_imp, imp_idx_nonzero, area_imp_final, xic_peak_rt_bounds, max_label, ratio_each_XIC_peak] = ...
    CIMPQuantCore.requantGroup(obj.m_ms12DatasetIO, group.rawName,...
    group.ratio(group.chargeGroupIdxs,:),...
    group.rts(group.chargeGroupIdxs,:),...
    group.intensity(group.chargeGroupIdxs,:),...
    group.lowMzBound, group.highMzBound, ...
    group.selectedCharge, group.impRtRanges,...
    obj.m_minMSMSnum);

if ~has_nonzero_imp
    return;
end
for i_imp = 1:length(imp_idx_nonzero)
    rt_peaks = struct(...
        'rt_start', xic_peak_rt_bounds(i_imp).rt_start, ...
        'rt_end', xic_peak_rt_bounds(i_imp).rt_end, ...
        'ratio', ratio_each_XIC_peak(i_imp), ...
        'check_label', max_label(i_imp));

    if isResearchRtMergeEnabled()
        rt_peaks = [rt_peaks, group.impRtRanges{imp_idx_nonzero(i_imp)}];
    end

    imp_records(end+1,1) = CIMPQuantRecord(...
        group.impNames{imp_idx_nonzero(i_imp)}, ...
        group.selectedCharge, ...
        group.rawName, ...
        mean([group.lowMzBound, group.highMzBound]), ...
        group.lowMzBound, ...
        group.highMzBound, ...
        area_imp_final(i_imp,1), ...
        rt_peaks); %#ok<AGROW>
end
end

function tf = isResearchRtMergeEnabled()
persistent cached_enabled has_cached_value
if isempty(has_cached_value)
    has_cached_value = false;
end

if ~has_cached_value
    env_val = strtrim(getenv('REQUANT_RT_PEAKS_MERGE_ON'));
    cached_enabled = any(strcmpi(env_val, {'1', 'true', 'yes', 'on'}));
    has_cached_value = true;
end

tf = cached_enabled;
end
