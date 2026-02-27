function idx_selected = select_by_alignment(imp_max_props, area_imp_by_peak, xic_peak_rt_bounds, rt_pred, options)
% Select peaks using alignment-aware scoring.
% Input:
%   imp_max_props (K x P double)
%       Max ratio contribution per IMP per peak
%   area_imp_by_peak (K x P double)
%       Area contribution per IMP per peak
%   xic_peak_rt_bounds (K x P struct)
%       Peak bounds per IMP
%   rt_pred (K x 1 double)
%       Predicted RT centers per IMP
%   options (struct, optional)
%       Selector options (rt_sigma)
% Output:
%   idx_selected (K x 1 double)
%       Selected peak index per IMP

if nargin < 5
    options = struct();
end

num_imp = size(imp_max_props, 1);
num_peaks = size(imp_max_props, 2);
if num_imp == 0 || num_peaks == 0
    idx_selected = zeros(num_imp, 1);
    return;
end

if numel(rt_pred) ~= num_imp
    error('select_by_alignment:InvalidRtPredSize', ...
        'rt_pred must have the same number of elements as imp_max_props rows.');
end

sigma = CStructOptionUtils.get(options, 'rt_sigma', 0.5);
weighted_props = imp_max_props;
for idx_imp = 1:num_imp
    rt_center = rt_pred(idx_imp);
    if isnan(rt_center)
        continue;
    end
    weights = CXICAlignScore.compute_rt_weights(rt_center, xic_peak_rt_bounds(idx_imp, :), sigma);
    weighted_props(idx_imp, :) = weighted_props(idx_imp, :) .* weights;
end

idx_selected = CXICPeakUtils.select_best_peak_per_imp(weighted_props, area_imp_by_peak);
end

