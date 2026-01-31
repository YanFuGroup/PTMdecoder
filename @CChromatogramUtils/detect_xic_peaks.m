function XIC_peaks = detect_xic_peaks(xic_rt, xic_intensity_smoothed, xic_intensity_raw, rt_sorted, alpha)
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
    % Output:
    %   XIC_peaks (1 x P struct)
    %       Struct array with fields: left_bound/right_bound (indices into xic_rt)
    
    % Extract the XIC peaks around the identified MSMS precursor
    % Extract from left to right
    idx_PSM = 1;
    XIC_peaks = struct('left_bound',{},'right_bound',{});
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
        right_bound = idx_first_rt;
        min_peak_inten = xic_intensity_smoothed(idx_first_rt);
        min_peak_iter = idx_first_rt;
        for iter_rt = idx_first_rt+1:length(xic_rt)
            % check whether reach the next rt of identified MS/MS precursor
            if idx_PSM<length(rt_sorted) && xic_rt(iter_rt)>rt_sorted(idx_PSM+1)
                right_bound = iter_rt;
                idx_PSM = idx_PSM + 1;
            end
            % update the local maximum
            if max_peak_inten < xic_intensity_smoothed(iter_rt)
                max_peak_inten = xic_intensity_smoothed(iter_rt);
                min_peak_inten = max_peak_inten;
                min_peak_iter = iter_rt;
            end
            % find the local minimum
            if xic_intensity_smoothed(iter_rt) < min_peak_inten
                min_peak_inten = xic_intensity_smoothed(iter_rt);
                min_peak_iter = iter_rt;
            end
            % find the right bound
            % two criteria: 1. the right is too low; 2. the local minimum is too low
            if xic_intensity_smoothed(iter_rt) < max_peak_inten*alpha
                right_bound = iter_rt - 1;
                break;
            elseif min_peak_inten < 0.5 * min(xic_intensity_smoothed(iter_rt),max_peak_inten)
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
        min_peak_inten = xic_intensity_smoothed(idx_first_rt);
        min_peak_iter = idx_first_rt;
        for iter_rt = idx_first_rt-1:-1:1
            % update the local maximum
            if max_peak_inten < xic_intensity_smoothed(iter_rt)
                max_peak_inten = xic_intensity_smoothed(iter_rt);
                % Update the min peak to ensure that the local minimum is after the maximum
                min_peak_inten = max_peak_inten;
                min_peak_iter = iter_rt;
            end
            % find the local minimum
            if xic_intensity_smoothed(iter_rt) < min_peak_inten
                min_peak_inten = xic_intensity_smoothed(iter_rt);
                min_peak_iter = iter_rt;
            end
            % find the left bound
            % two criteria: 1. the left is too low; 2. the minimum is too low
            if xic_intensity_smoothed(iter_rt) < max_peak_inten*alpha
                XIC_peaks(i_Xp).left_bound = iter_rt + 1;
                break;
            elseif min_peak_inten < 0.5 * min(xic_intensity_smoothed(iter_rt),max_peak_inten)
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
        if sum(xic_intensity_raw(XIC_peaks(i_Xp).left_bound:XIC_peaks(i_Xp).right_bound)~=0) < 5
            XIC_peaks(i_Xp) = [];
        end
    end
end
