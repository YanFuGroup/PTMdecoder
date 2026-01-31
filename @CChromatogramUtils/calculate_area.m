function area = calculate_area(xic_rt, xic_intensity_full, idx_start, idx_end)
% calculate_area
% Calculates the area under the curve (AUC) for a peak using closed integration.
% Internally uses get_closed_peak_data to ensure the peak is "closed" (starts and ends at 0).
%
% Input:
%   xic_rt (N x 1 double) minutes
%       Vector of retention times
%   xic_intensity_full (N x 1 double) intensity
%       Vector of intensities corresponding to xic_rt
%   idx_start (1 x 1 double/int)
%       Start index of the peak (internal)
%   idx_end (1 x 1 double/int)
%       End index of the peak (internal)
%
% Output:
%   area (1 x 1 double) intensity * seconds
%       Calculated area (trapz * 60; assumes RT in minutes)
%
% Attention:
%   xic_rt in in minutes, area output in intensity * seconds.
%   This function assumes that the input indices correspond to a peak
%   that has been properly identified and that the intensity vector is
%   aligned with the xic_rt.

    % Get closed peak data (handling boundary expansion and zero-padding)
    [rec_rt, rec_inten] = CChromatogramUtils.get_closed_peak_data(...
        xic_rt, xic_intensity_full, idx_start, idx_end);
    
    % Calculate area using trapezoidal integration
    % Multiplied by 60 to convert area units (assuming RT in minutes -> output in intensity * seconds)
    area = trapz(rec_rt, rec_inten) * 60;
end
