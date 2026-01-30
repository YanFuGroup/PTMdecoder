function peak_ranges = map_rt_to_indices(rt_grid, final_XIC_peak, is_skip_vec, rt_error_tol)
% map_rt_to_indices
% Maps retention time bounds to indices in the retention time grid.
% Each IMP's RT bounds are converted to corresponding index bounds in rt_grid.
%
% Input:
%   rt_grid:        Vector of retention times (the grid)
%   final_XIC_peak: Structure array with .left_bound and .right_bound (RT values)
%   is_skip_vec:    Boolean vector indicating which IMP to skip
%   rt_error_tol:   Tolerance for finding the specific retention time
%
% Output:
%   peak_ranges:    Structure array with .left_bound and .right_bound (Index values)

    num_imp = length(final_XIC_peak);
    peak_ranges = repmat(struct('left_bound',0,'right_bound',0), num_imp, 1);
    
    for idx_imp = 1:num_imp
        if is_skip_vec(idx_imp), continue; end
        
        [diff_l, peak_ranges(idx_imp).left_bound] = min(abs(rt_grid - final_XIC_peak(idx_imp).left_bound));
        if diff_l > rt_error_tol
            error(['Cannot find the spectra on the specified retention time: ', ...
                num2str(final_XIC_peak(idx_imp).left_bound)]);
        end
        
        [diff_r, peak_ranges(idx_imp).right_bound] = min(abs(rt_grid - final_XIC_peak(idx_imp).right_bound));
        if diff_r > rt_error_tol
            error(['Cannot find the spectra on the specified retention time: ', ...
                num2str(final_XIC_peak(idx_imp).right_bound)]);
        end
    end
end
