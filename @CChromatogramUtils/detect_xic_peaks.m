function XIC_peaks = detect_xic_peaks(rt_grid, smoothed_intensity, raw_intensity, sort_rts, alpha)
    % Detect XIC peaks using smoothed intensity and identified MSMS events
    % Inputs:
    %   rt_grid: Retention time grid
    %   smoothed_intensity: Smoothed XIC
    %   raw_intensity: Raw XIC (for filtering small peaks)
    %   sort_rts: Sorted retention times of PSMs
    %   alpha: Threshold parameter for peak boundary detection
    
    % Extract the XIC peaks around the identified MSMS precursor
    % Extract from left to right
    idx_PSM = 1;
    XIC_peaks = struct('left_bound',{},'right_bound',{});
    i_Xp = 1;
    while idx_PSM <= length(sort_rts)
        % index of rt for first identified MS/MS in this peak
        [~, idx_first_rt] = min(abs(rt_grid-sort_rts(idx_PSM)));
        max_peak_inten = smoothed_intensity(idx_first_rt);
        if max_peak_inten == 0
            % An identification with intensity of zero means mis-identification
            idx_PSM = idx_PSM + 1;
            continue;
        end

        % look after the right boundary
        right_bound = idx_first_rt;
        min_peak_inten = smoothed_intensity(idx_first_rt);
        min_peak_iter = idx_first_rt;
        for iter_rt = idx_first_rt+1:length(rt_grid)
            % check whether reach the next rt of identified MS/MS precursor
            if idx_PSM<length(sort_rts) && rt_grid(iter_rt)>sort_rts(idx_PSM+1)
                right_bound = iter_rt;
                idx_PSM = idx_PSM + 1;
            end
            % update the local maximum
            if max_peak_inten < smoothed_intensity(iter_rt)
                max_peak_inten = smoothed_intensity(iter_rt);
                min_peak_inten = max_peak_inten;
                min_peak_iter = iter_rt;
            end
            % find the local minimum
            if smoothed_intensity(iter_rt) < min_peak_inten
                min_peak_inten = smoothed_intensity(iter_rt);
                min_peak_iter = iter_rt;
            end
            % find the right bound
            % two criteria: 1. the right is too low; 2. the local minimum is too low
            if smoothed_intensity(iter_rt) < max_peak_inten*alpha
                right_bound = iter_rt - 1;
                break;
            elseif min_peak_inten < 0.5 * min(smoothed_intensity(iter_rt),max_peak_inten)
                right_bound = min_peak_iter;
                break;
            end
        end
        
        % Check index bounds
        if i_Xp > length(XIC_peaks)
            XIC_peaks(i_Xp).left_bound = 0; % Initialize
        end 
        XIC_peaks(i_Xp).right_bound = right_bound;

        % look after the left boundary
        XIC_peaks(i_Xp).left_bound = idx_first_rt;
        min_peak_inten = smoothed_intensity(idx_first_rt);
        min_peak_iter = idx_first_rt;
        for iter_rt = idx_first_rt-1:-1:1
            % update the local maximum
            if max_peak_inten < smoothed_intensity(iter_rt)
                max_peak_inten = smoothed_intensity(iter_rt);
                % Update the min peak to ensure that the local minimum is after the maximum
                min_peak_inten = max_peak_inten;
                min_peak_iter = iter_rt;
            end
            % find the local minimum
            if smoothed_intensity(iter_rt) < min_peak_inten
                min_peak_inten = smoothed_intensity(iter_rt);
                min_peak_iter = iter_rt;
            end
            % find the left bound
            % two criteria: 1. the left is too low; 2. the minimum is too low
            if smoothed_intensity(iter_rt) < max_peak_inten*alpha
                XIC_peaks(i_Xp).left_bound = iter_rt + 1;
                break;
            elseif min_peak_inten < 0.5 * min(smoothed_intensity(iter_rt),max_peak_inten)
                XIC_peaks(i_Xp).left_bound = min_peak_iter;
                break;
            end
        end

        % prepare for next peak extraction
        i_Xp = i_Xp + 1;
        idx_PSM = idx_PSM + 1;
    end

    % remove the XIC peaks only with less than 5 point
    for i_Xp = length(XIC_peaks):-1:1
        if sum(raw_intensity(XIC_peaks(i_Xp).left_bound:XIC_peaks(i_Xp).right_bound)~=0) < 5
            XIC_peaks(i_Xp) = [];
        end
    end
end
