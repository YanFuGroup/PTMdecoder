function weights = compute_rt_weights(rt_pred, xic_peak_rt_bounds, sigma)
% Compute RT weights for peak selection.
% Input:
%   rt_pred (1 x 1 double)
%       Predicted RT center (minutes)
%   xic_peak_rt_bounds (1 x P struct)
%       Peak bounds with rt_start/rt_end
%   sigma (1 x 1 double, optional)
%       Decay factor for RT residuals
% Output:
%   weights (1 x P double)
%       Weight per peak based on RT residual

if nargin < 3 || isempty(sigma)
    sigma = 0.5;
end

num_peaks = numel(xic_peak_rt_bounds);
weights = zeros(1, num_peaks);
for idx = 1:num_peaks
    center_rt = mean([xic_peak_rt_bounds(idx).rt_start, xic_peak_rt_bounds(idx).rt_end]);
    residual = abs(center_rt - rt_pred);
    weights(idx) = exp(-residual / sigma);
end
end
