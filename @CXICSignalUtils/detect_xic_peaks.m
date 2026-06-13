function [xic_peak_idx_bounds, diagnostics] = detect_xic_peaks(xic_rt, xic_intensity_smoothed, xic_intensity_raw, rt_sorted, alpha, min_nonzero_points)
    % Detect XIC peaks using smoothed intensity and identified MSMS events
    % Inputs:
    %   xic_rt (N x 1 double) minutes
    %       Retention time grid
    %   xic_intensity_smoothed (N x 1 double) intensity
    %       Smoothed XIC
    %   xic_intensity_raw (N x 1 double) intensity
    %       Raw XIC (for filtering small peaks)
    %   rt_sorted (M x 1 double) minutes
    %       Sorted retention times of PSMs
    %   alpha (1 x 1 double)
    %       Threshold parameter for peak boundary detection
    %   min_nonzero_points (1 x 1 double/int, optional)
    %       minimum non-zero raw MS1 points required within a candidate peak
    % Output:
    %   xic_peak_idx_bounds (1 x P struct)
    %       Struct array with fields: idx_start/idx_end (indices into xic_rt)
    
    if nargin < 6 || isempty(min_nonzero_points)
        min_nonzero_points = 5;
    end
    diagnostics = struct( ...
        'candidate_peak_count', 0, ...
        'filtered_sparse_peak_count', 0, ...
        'candidate_nonzero_points', zeros(0, 1), ...
        'min_nonzero_points', min_nonzero_points);

    % Extract the XIC peaks around the identified MSMS precursor
    % Extract from left to right
    idx_PSM = 1;
    xic_peak_idx_bounds = struct('idx_start',{},'idx_end',{});
    i_Xp = 1;
    while idx_PSM <= length(rt_sorted)
        % index of rt for first identified MS/MS in this peak
        [~, idx_first_rt] = min(abs(xic_rt-rt_sorted(idx_PSM)));
        max_peak_inten = xic_intensity_smoothed(idx_first_rt);
        if max_peak_inten == 0
            % An identification with intensity of zero means mis-identification
            idx_PSM = idx_PSM + 1;
            continue;
        end

        % look after the right boundary
        idx_end = idx_first_rt;
        min_peak_inten = xic_intensity_smoothed(idx_first_rt);
        min_peak_iter = idx_first_rt;
        for idx_xic_rt = idx_first_rt+1:length(xic_rt)
            % check whether reach the next rt of identified MS/MS precursor
            if idx_PSM<length(rt_sorted) && xic_rt(idx_xic_rt)>rt_sorted(idx_PSM+1)
                idx_end = idx_xic_rt;
                idx_PSM = idx_PSM + 1;
            end
            % update the local maximum
            if max_peak_inten < xic_intensity_smoothed(idx_xic_rt)
                max_peak_inten = xic_intensity_smoothed(idx_xic_rt);
                min_peak_inten = max_peak_inten;
                min_peak_iter = idx_xic_rt;
            end
            % find the local minimum
            if xic_intensity_smoothed(idx_xic_rt) < min_peak_inten
                min_peak_inten = xic_intensity_smoothed(idx_xic_rt);
                min_peak_iter = idx_xic_rt;
            end
            % find the right bound
            % two criteria: 1. the right is too low; 2. the local minimum is too low
            if xic_intensity_smoothed(idx_xic_rt) < max_peak_inten*alpha
                idx_end = idx_xic_rt - 1;
                break;
            elseif min_peak_inten < 0.5 * min(xic_intensity_smoothed(idx_xic_rt),max_peak_inten)
                idx_end = min_peak_iter;
                break;
            end
        end
        
        % Check index bounds
        if i_Xp > length(xic_peak_idx_bounds)
            xic_peak_idx_bounds(i_Xp).idx_start = 0; % Initialize
        end 
        xic_peak_idx_bounds(i_Xp).idx_end = idx_end;

        % look after the left boundary
        xic_peak_idx_bounds(i_Xp).idx_start = idx_first_rt;
        min_peak_inten = xic_intensity_smoothed(idx_first_rt);
        min_peak_iter = idx_first_rt;
        for idx_xic_rt = idx_first_rt-1:-1:1
            % update the local maximum
            if max_peak_inten < xic_intensity_smoothed(idx_xic_rt)
                max_peak_inten = xic_intensity_smoothed(idx_xic_rt);
                % Update the min peak to ensure that the local minimum is after the maximum
                min_peak_inten = max_peak_inten;
                min_peak_iter = idx_xic_rt;
            end
            % find the local minimum
            if xic_intensity_smoothed(idx_xic_rt) < min_peak_inten
                min_peak_inten = xic_intensity_smoothed(idx_xic_rt);
                min_peak_iter = idx_xic_rt;
            end
            % find the left bound
            % two criteria: 1. the left is too low; 2. the minimum is too low
            if xic_intensity_smoothed(idx_xic_rt) < max_peak_inten*alpha
                xic_peak_idx_bounds(i_Xp).idx_start = idx_xic_rt + 1;
                break;
            elseif min_peak_inten < 0.5 * min(xic_intensity_smoothed(idx_xic_rt),max_peak_inten)
                xic_peak_idx_bounds(i_Xp).idx_start = min_peak_iter;
                break;
            end
        end

        % prepare for next peak extraction
        i_Xp = i_Xp + 1;
        idx_PSM = idx_PSM + 1;
    end

    diagnostics.candidate_peak_count = length(xic_peak_idx_bounds);
    diagnostics.candidate_nonzero_points = zeros(diagnostics.candidate_peak_count, 1);

    % Remove XIC peaks with too few non-zero raw MS1 points.
    for i_Xp = length(xic_peak_idx_bounds):-1:1
        nonzero_points = sum(xic_intensity_raw( ...
            xic_peak_idx_bounds(i_Xp).idx_start:xic_peak_idx_bounds(i_Xp).idx_end) ~= 0);
        diagnostics.candidate_nonzero_points(i_Xp) = nonzero_points;
        if nonzero_points < min_nonzero_points
            diagnostics.filtered_sparse_peak_count = diagnostics.filtered_sparse_peak_count + 1;
            xic_peak_idx_bounds(i_Xp) = [];
        end
    end
end
