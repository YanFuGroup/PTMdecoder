function [final_XIC_peak_for_IMP, max_label, is_skip_vec, peak_ranges] = ...
    prepare_peak_ranges_from_imp_rt_range(rt_grid, current_iso_rt_range, rt_error_tol)
% Prepare peak ranges based on input iso RT ranges.
% input:
%   rt_grid (N x 1 double) minutes
%       retention time grid
%   current_iso_rt_range (K x 1 cell)
%       RT ranges for each IMP (each cell: [] or [start, end] in minutes)
%   rt_error_tol (1 x 1 double) minutes
%       RT tolerance in minutes
% output:
%   final_XIC_peak_for_IMP (K x 1 struct)
%       RT bounds for each IMP peak; fields: left_bound/right_bound (minutes)
%   max_label (K x 1 double)
%       max check label for each IMP
%   is_skip_vec (K x 1 logical)
%       vector indicating IMPs to skip
%   peak_ranges (K x 1 struct)
%       index bounds for each IMP peak; fields: left_bound/right_bound (indices into rt_grid)

is_skip_vec = cellfun(@isempty, current_iso_rt_range);

[final_XIC_peak_for_IMP, max_label, is_skip_vec] = ...
    CChromatogramUtils.parse_imp_rt_ranges(current_iso_rt_range, is_skip_vec);

peak_ranges = CChromatogramUtils.map_rt_to_indices(...
    rt_grid, final_XIC_peak_for_IMP, is_skip_vec, rt_error_tol);
end
