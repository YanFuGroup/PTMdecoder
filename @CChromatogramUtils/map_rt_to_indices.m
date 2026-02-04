function xic_peak_idx_bounds = map_rt_to_indices(xic_rt, final_xic_peak_rt_bounds, is_skip_vec, rt_error_tol)
% map_rt_to_indices
% Maps retention time bounds to indices in the retention time grid.
% Each IMP's RT bounds are converted to corresponding index bounds in xic_rt.
%
% Input:
%   xic_rt (N x 1 double) minutes
%       Vector of retention times (the grid)
%   final_xic_peak_rt_bounds (K x 1 struct)
%       Structure array with .rt_start and .rt_end (RT values, minutes)
%   is_skip_vec (K x 1 logical)
%       Boolean vector indicating which IMP to skip
%   rt_error_tol (1 x 1 double) minutes
%       Tolerance for finding the specific retention time
%
% Output:
%   xic_peak_idx_bounds (K x 1 struct)
%       Structure array with .idx_start and .idx_end (index values into xic_rt)

    num_imp = length(final_xic_peak_rt_bounds);
    xic_peak_idx_bounds = repmat(struct('idx_start',0,'idx_end',0), num_imp, 1);
    
    for idx_imp = 1:num_imp
        if is_skip_vec(idx_imp), continue; end
        
        [diff_l, xic_peak_idx_bounds(idx_imp).idx_start] = min(abs(xic_rt - final_xic_peak_rt_bounds(idx_imp).rt_start));
        if diff_l > rt_error_tol
            error(['Cannot find the spectra on the specified retention time: ', ...
            num2str(final_xic_peak_rt_bounds(idx_imp).rt_start)]);
        end
        
        [diff_r, xic_peak_idx_bounds(idx_imp).idx_end] = min(abs(xic_rt - final_xic_peak_rt_bounds(idx_imp).rt_end));
        if diff_r > rt_error_tol
            error(['Cannot find the spectra on the specified retention time: ', ...
            num2str(final_xic_peak_rt_bounds(idx_imp).rt_end)]);
        end
    end
end
