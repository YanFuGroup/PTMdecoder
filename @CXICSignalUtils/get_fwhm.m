function fwhm = get_fwhm(peak_rts, peak_intens)
    % Get the full width at half maxima (FWHM)
    % Input:
    %   peak_rts (N x 1 double) minutes
    %       Retention times of peaks
    %   peak_intens (N x 1 double) intensity
    %       Intensities of peaks
    % Output:
    %   fwhm (1 x 1 double) minutes
    %       Full width at half maxima
    
    % Default, used for zero-intensity peak or null peak
    fwhm = 0;
    
    % Check whether the lengths of peak_rts and peak_intens are equal
    if length(peak_rts) ~= length(peak_intens)
        error('The length of peak_rts (%d) and peak_intens (%d) is not equal',...
            length(peak_rts), length(peak_intens));
    end
    
    % find two ranges of half maxima
    half_maxima = max(peak_intens)/2;
    if half_maxima==0
        return;
    end
    % the two point around the left half maxima, find left half bound
    left_right = find(peak_intens>half_maxima, 1);
    left_left = left_right - 1;
    if left_left > 0
        left_half_point = peak_rts(left_left)+(half_maxima-peak_intens(left_left))*...
            (peak_rts(left_right)-peak_rts(left_left))/(peak_intens(left_right)-peak_intens(left_left));
    else
        % when the left half maxima is smaller than the most left point in
        %   this peak, use the most left point to measure the peak width
        left_half_point = peak_rts(1);
    end
    % the two point around the right half maxima, find right half bound
    right_left = find(peak_intens>half_maxima, 1, 'last');
    right_right = right_left + 1;
    if right_right < length(peak_rts)
        right_half_point = peak_rts(right_right)-(half_maxima-peak_intens(right_right))*...
            (peak_rts(right_right)-peak_rts(right_left))/(peak_intens(right_left)-peak_intens(right_right));
    else
        % when the left half maxima is smaller than the most left point in
        %   this peak, use the most left point to measure the peak width
        right_half_point = peak_rts(end);
    end
    fwhm = right_half_point - left_half_point;
end
