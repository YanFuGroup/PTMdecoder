function esti_ratio = calculate_kernel_ratio(rt_grid, sort_rts, sort_ratioMatrix, peak_ranges, is_broadcast)
    % Calculate estimated ratio using Gaussian kernel
    % Inputs:
    %   rt_grid: RT grid (vector)
    %   sort_rts: Sorted PSM RTs
    %   sort_ratioMatrix: quantification matrix
    %   peak_ranges: struct array with .left_bound and .right_bound (indices)
    %   is_broadcast: boolean, true if every peak_range applies to all IMPs (Quant), 
    %              false if each peak corresponds to each IMP (Requant)
    
    num_imp = size(sort_ratioMatrix, 2);
    esti_ratio = zeros(length(rt_grid), num_imp);
    
    % Gaussian kernel function
    kernal_func = @(u) (1/sqrt(2*pi))*exp(-0.5*u.^2);
    
    % Determine loop range and validation
    if is_broadcast
        % Quant mode: Iterate all detected peaks
        n_loops = length(peak_ranges);
    else
        % Requant mode: One-to-one mapping
        % Strict Validation: peak_ranges count MUST match IMP count
        if length(peak_ranges) ~= num_imp
             error('CChromatogramUtils:InvalidInput', ...
                   'In Requant mode (is_broadcast=false), peak_ranges length (%d) must match number of IMPs (%d).', ...
                   length(peak_ranges), num_imp);
        end
        n_loops = num_imp;
    end

    % Unified calculation loop
    for i = 1:n_loops
        % 1. Prepare Data
        % Range is from the i-th peak range structure (guaranteed valid by n_loops definition)
        left = peak_ranges(i).left_bound;
        right = peak_ranges(i).right_bound;
        
        range_indices = left:right;
        
        % Common validity checks
        if isempty(range_indices) || isempty(left) || isempty(right) 
            continue;
        end
        if ~is_broadcast && isequal(left, right) % Specific check often seen in requant
            continue;
        end
        
        % Collect rts within current peak
        eps_rt_print = 1e-4;
        idxs_rt_mask = sort_rts>=rt_grid(left)-eps_rt_print & sort_rts<=rt_grid(right)+eps_rt_print;
        rts_current = sort_rts(idxs_rt_mask);
        
        if isempty(rts_current), continue; end
        
        % 2. Core Calculation (Bandwidth & Weights)
        bandwidth_val = (4/(3*size(rts_current,1)))^0.2*std(rts_current);
        weights = zeros(length(range_indices), length(rts_current));
        
        for idx_PSM = 1:length(rts_current)
            if bandwidth_val == 0, break; end
            weights(:, idx_PSM) = kernal_func((rt_grid(range_indices) - rts_current(idx_PSM))/bandwidth_val);
        end
        
        % Handle zero bandwidth or zero weights
        if bandwidth_val == 0 || any(all(weights<1e-15/length(sort_rts), 2))
            bandwidth_val = min(rt_grid(right)-rt_grid(left), 1);
            for idx_PSM = 1:length(rts_current)
                weights(:, idx_PSM) = kernal_func((rt_grid(range_indices) - rts_current(idx_PSM))/bandwidth_val);
            end
        end
        
        % 3. Apply Result (Policy Decision)
        if is_broadcast
            target_indices = 1:num_imp; % Broadcast result to ALL IMPs
        else
            target_indices = i;         % Apply result ONLY to the current i-th IMP
        end
        
        weight_sum = sum(weights, 2) + eps;
        for idx_iso = target_indices
            esti_ratio(range_indices, idx_iso) = esti_ratio(range_indices, idx_iso) + ...
                (weights * sort_ratioMatrix(idxs_rt_mask, idx_iso)) ./ weight_sum;
        end
    end
    
    % Normalize
    esti_ratio = esti_ratio ./ (sum(esti_ratio, 2) + eps);
end
