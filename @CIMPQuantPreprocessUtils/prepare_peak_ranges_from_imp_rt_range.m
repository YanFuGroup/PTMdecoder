function [xic_peak_rt_bounds, max_label, is_skip_vec, xic_peak_idx_bounds] = ...
    prepare_peak_ranges_from_imp_rt_range(xic_rt, current_imp_rt_range, rt_error_tol)
% Prepare peak ranges based on input iso RT ranges.
% input:
%   xic_rt (N x 1 double) minutes
%       retention time grid
%   current_imp_rt_range (K x 1 cell)
%       RT ranges for each IMP (each cell: [] or [start, end] in minutes)
%   rt_error_tol (1 x 1 double) minutes
%       RT tolerance in minutes
% output:
%   xic_peak_rt_bounds (K x 1 struct)
%       RT bounds for each IMP peak; fields: rt_start/rt_end (minutes)
%   max_label (K x 1 double)
%       max check label for each IMP
%   is_skip_vec (K x 1 logical)
%       vector indicating IMPs to skip
%   xic_peak_idx_bounds (K x 1 struct)
%       index bounds for each IMP peak; fields: idx_start/idx_end (indices into xic_rt)

is_skip_vec = cellfun(@isempty, current_imp_rt_range);

[xic_peak_rt_bounds, max_label, is_skip_vec] = ...
    CChromatogramUtils.parse_imp_rt_ranges(current_imp_rt_range, is_skip_vec);

xic_peak_idx_bounds = CChromatogramUtils.map_rt_to_indices(...
    xic_rt, xic_peak_rt_bounds, is_skip_vec, rt_error_tol);
end
