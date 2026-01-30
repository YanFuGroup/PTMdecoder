function [final_XIC_peak_for_IMP, max_label, is_skip_vec, peak_ranges] = ...
    prepare_peak_ranges_from_iso_rt_range(rt_grid, current_iso_rt_range, rt_error_tol)
% Prepare peak ranges based on input iso RT ranges.
% input:
%   rt_grid
%       retention time grid
%   current_iso_rt_range
%       retention times of current IMPs
%   rt_error_tol
%       RT tolerance in minutes
% output:
%   final_XIC_peak_for_IMP
%       RT bounds for each IMP peak
%   max_label
%       max check label for each IMP
%   is_skip_vec
%       vector indicating IMPs to skip
%   peak_ranges
%       index bounds for each IMP peak

is_skip_vec = cellfun(@isempty, current_iso_rt_range);

[final_XIC_peak_for_IMP, max_label, is_skip_vec] = ...
    CChromatogramUtils.parse_imp_rt_ranges(current_iso_rt_range, is_skip_vec);

peak_ranges = CChromatogramUtils.map_rt_to_indices(...
    rt_grid, final_XIC_peak_for_IMP, is_skip_vec, rt_error_tol);
end
