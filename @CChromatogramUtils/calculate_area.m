function area = calculate_area(rt_grid, intensity_full, idx_start, idx_end)
% calculate_area
% Calculates the area under the curve (AUC) for a peak using closed integration.
% Internally uses get_closed_peak_data to ensure the peak is "closed" (starts and ends at 0).
%
% Input:
%   rt_grid:        Vector of retention times
%   intensity_full: Vector of intensities corresponding to rt_grid
%   idx_start:      Start index of the peak (internal)
%   idx_end:        End index of the peak (internal)
%
% Output:
%   area:           Calculated area (using trapz * 60 for minutes to seconds conversion if applicable) 
%                   Note: The * 60 factor assumes RT is in minutes.
% 
% Attention:
%   rt_grid in in minutes, area output in intensity * seconds.
%   This function assumes that the input indices correspond to a peak
%   that has been properly identified and that the intensity vector is
%   aligned with the rt_grid.

    % Get closed peak data (handling boundary expansion and zero-padding)
    [rec_rt, rec_inten] = CChromatogramUtils.get_closed_peak_data(...
        rt_grid, intensity_full, idx_start, idx_end);
    
    % Calculate area using trapezoidal integration
    % Multiplied by 60 to convert area units (assuming RT in minutes -> output in intensity * seconds)
    area = trapz(rec_rt, rec_inten) * 60;
end
