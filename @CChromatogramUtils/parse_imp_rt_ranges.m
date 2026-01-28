function [final_XIC_peak, max_label, is_skip_vec] = parse_imp_rt_ranges(imp_rt_range, is_skip_vec)
% parse_imp_rt_ranges
% Parses the retention time ranges structure to find the best peak for each IMP.
%
% Input:
%   imp_rt_range: Cell array of structures containing peak info (check_label, rt_start, rt_end)
%   is_skip_vec:  Boolean vector indicating which IMP to skip
%
% Output:
%   final_XIC_peak: Structure array with .left_bound and .right_bound (RT values)
%   max_label:      Vector of max check labels found for each IMP
%   is_skip_vec:    Updated boolean vector (skips if max_label is 0)

    num_imp = length(imp_rt_range);
    final_XIC_peak = repmat(struct('left_bound',0,'right_bound',0), num_imp, 1);
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
        
        final_XIC_peak(idx_imp).left_bound = curr_range(idx_max).rt_start;
        final_XIC_peak(idx_imp).right_bound = curr_range(idx_max).rt_end;
    end
end
