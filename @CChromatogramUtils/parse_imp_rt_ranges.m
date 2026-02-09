function [xic_peak_rt_bounds, max_label, is_skip_vec] = parse_imp_rt_ranges(imp_rt_range, is_skip_vec)
% parse_imp_rt_ranges
% Parses the retention time ranges structure to find the best peak for each IMP.
%
% Input:
%   imp_rt_range (K x 1 cell)
%       Cell array of structures containing peak info (check_label, rt_start, rt_end)
%   is_skip_vec (K x 1 logical)
%       Boolean vector indicating which IMP to skip
%
% Output:
%   xic_peak_rt_bounds (K x 1 struct)
%       Structure array with .rt_start and .rt_end (RT values in minutes)
%   max_label (K x 1 double)
%       Vector of max check labels found for each IMP
%   is_skip_vec (K x 1 logical)
%       Updated boolean vector (skips if max_label is 0)

    num_imp = length(imp_rt_range);
    xic_peak_rt_bounds = repmat(struct('rt_start',0,'rt_end',0), num_imp, 1);
    max_label = zeros(num_imp,1);

    for idx_imp = 1:num_imp
        % Check if need to consider this IMP
        if is_skip_vec(idx_imp)
            continue;
        end

        curr_range = imp_rt_range{idx_imp};
        if isempty(curr_range)
             is_skip_vec(idx_imp) = true;
             continue;
        end

        % Record all of the check labels
        check_labels = [curr_range.check_label];
        
        % Find the peak with max check label (the first of max peaks)
        [max_label(idx_imp), idx_max] = max(check_labels);
        
        if max_label(idx_imp) == 0
            % If all check labels are zero, skip this IMP later on
            is_skip_vec(idx_imp) = true;
            continue;
        end
        
        xic_peak_rt_bounds(idx_imp).rt_start = curr_range(idx_max).rt_start;
        xic_peak_rt_bounds(idx_imp).rt_end = curr_range(idx_max).rt_end;
    end
end
