function [imp_max_props, peak_fwhms, area_imp_by_peak, xic_peak_rt_bounds] = compute_peak_features(xic_rt, xic_intensity_smoothed, xic_ratio_estimated, xic_peak_idx_bounds)
% compute_peak_features
% Compute IMP-wise peak features across all candidate XIC peaks.
%
% Inputs:
%   xic_rt (N x 1 double) minutes
%       RT grid vector
%   xic_intensity_smoothed (N x 1 double) intensity
%       Smoothed XIC intensity (aligned to xic_rt)
%   xic_ratio_estimated (N x K double)
%       Estimated ratio matrix for K IMPs at each RT; rows sum to ~1 within peaks
%   xic_peak_idx_bounds (1 x P struct)
%       Struct array with fields: idx_start/idx_end (indices into xic_rt)
%
% Outputs:
%   imp_max_props (K x P double)
%       Max ratio contribution per IMP per peak
%   peak_fwhms (K x P double) minutes
%       FWHM per IMP per peak
%   area_imp_by_peak (K x P double) area
%       Area contribution per IMP per peak
%   xic_peak_rt_bounds (K x P struct)
%       RT bounds per IMP per peak, fields: .rt_start/.rt_end (minutes)

num_imp = size(xic_ratio_estimated, 2);
num_peaks = length(xic_peak_idx_bounds);

intensityMatrix = xic_ratio_estimated .* xic_intensity_smoothed;
area_imp_by_peak = zeros(num_imp, num_peaks);
peak_fwhms = zeros(num_imp, num_peaks);
imp_max_props = zeros(num_imp, num_peaks);

single_rt_bounds = repmat(struct('rt_start', 0, 'rt_end', 0), 1, num_peaks);

for i_Xp = 1:num_peaks
    curr_idx_start = xic_peak_idx_bounds(i_Xp).idx_start;
    curr_idx_end = xic_peak_idx_bounds(i_Xp).idx_end;
    
    % 1. XIC-dependent properties (Independent of IMP)
    peak_rts = xic_rt(curr_idx_start:curr_idx_end);
    
    single_rt_bounds(i_Xp).rt_start = xic_rt(curr_idx_start);
    single_rt_bounds(i_Xp).rt_end = xic_rt(curr_idx_end);
    
    % 2. IMP-dependent properties
    ratio_slice = xic_ratio_estimated(curr_idx_start:curr_idx_end, :);
    imp_max_props(:, i_Xp) = max(ratio_slice, [], 1)';
    
    % Calculate area and FWHM per IMP
    for idx_imp = 1:num_imp
        peak_fwhms(idx_imp, i_Xp) = CXICSignalUtils.get_fwhm( ...
            peak_rts, intensityMatrix(curr_idx_start:curr_idx_end, idx_imp));
        
        area_imp_by_peak(idx_imp, i_Xp) = CXICSignalUtils.calculate_area(...
            xic_rt, intensityMatrix(:, idx_imp), curr_idx_start, curr_idx_end);
    end
end

xic_peak_rt_bounds = repmat(single_rt_bounds, num_imp, 1);
end
