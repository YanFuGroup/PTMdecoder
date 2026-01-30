classdef CQuantIMPGroupUtils
    % CQuantIMPGroupUtils
    % Utilities for IMP group quantification workflows.
    %
    % This class hosts mid-level quantification logic that composes
    % chromatogram utilities into higher-level IMP-group behaviors.
    methods (Static)
        esti_ratio = filter_and_normalize_peak_ratios(rt_grid, smoothed_intensity, esti_ratio, XIC_peaks, resFilterThres)
    end
end
