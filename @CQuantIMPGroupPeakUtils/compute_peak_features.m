function [imp_max_props, peak_fwhms, area_each_XIC_peak, rt_bound] = compute_peak_features(rt_grid, smoothed_intensity, esti_ratio, XIC_peaks)
% compute_peak_features
% Compute IMP-wise peak features across all candidate XIC peaks.
%
% Inputs:
%   rt_grid (N x 1 double) minutes
%       RT grid vector
%   smoothed_intensity (N x 1 double) intensity
%       Smoothed XIC intensity (aligned to rt_grid)
%   esti_ratio (N x K double)
%       Estimated ratio matrix for K IMPs at each RT; rows sum to ~1 within peaks
%   XIC_peaks (1 x P struct)
%       Struct array with fields: left_bound/right_bound (indices into rt_grid)
%
% Outputs:
%   imp_max_props (K x P double)
%       Max ratio contribution per IMP per peak
%   peak_fwhms (K x P double) minutes
%       FWHM per IMP per peak
%   area_each_XIC_peak (K x P double) area
%       Area contribution per IMP per peak
%   rt_bound (K x P struct)
%       RT bounds per IMP per peak, fields: .start/.end (minutes)

num_imp = size(esti_ratio, 2);
num_peaks = length(XIC_peaks);

intensityMatrix = esti_ratio .* smoothed_intensity;
area_each_XIC_peak = zeros(num_imp, num_peaks);
peak_fwhms = zeros(num_imp, num_peaks);
imp_max_props = zeros(num_imp, num_peaks);

single_rt_bounds = repmat(struct('start', 0, 'end', 0), 1, num_peaks);

for i_Xp = 1:num_peaks
    curr_start = XIC_peaks(i_Xp).left_bound;
    curr_end = XIC_peaks(i_Xp).right_bound;
    
    % 1. XIC-dependent properties (Independent of IMP)
    peak_rts = rt_grid(curr_start:curr_end);
    
    single_rt_bounds(i_Xp).start = rt_grid(curr_start);
    single_rt_bounds(i_Xp).end = rt_grid(curr_end);
    
    % 2. IMP-dependent properties
    ratio_slice = esti_ratio(curr_start:curr_end, :);
    imp_max_props(:, i_Xp) = max(ratio_slice, [], 1)';
    
    % Calculate area and FWHM per IMP
    for idx_imp = 1:num_imp
        peak_fwhms(idx_imp, i_Xp) = CChromatogramUtils.get_fwhm( ...
            peak_rts, intensityMatrix(curr_start:curr_end, idx_imp));
        
        area_each_XIC_peak(idx_imp, i_Xp) = CChromatogramUtils.calculate_area(...
            rt_grid, intensityMatrix(:, idx_imp), curr_start, curr_end);
    end
end

rt_bound = repmat(single_rt_bounds, num_imp, 1);
end
