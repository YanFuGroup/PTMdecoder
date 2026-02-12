function model = fit(ref_rts, target_rts, options)
% Fit a global + local RT alignment transform.
% Input:
%   ref_rts (N x 1 double)
%       Reference RTs (minutes)
%   target_rts (N x 1 double)
%       Target RTs (minutes)
%   options (struct, optional)
%       Transform options (num_bins, min_per_bin)
%       num_bins: number of RT bins for local residual correction
%       min_per_bin: minimum points required per bin to estimate offset
%       outlier_k: k for outlier removal (0 to disable)
%       outlier_method: method for outlier detection ('mad' or 'std')
% Output:
%   model (struct)
%       Transform model with fields slope, intercept, bin_centers, bin_offsets, has_model

if nargin < 3
    options = struct();
end

model = struct('slope', 1, 'intercept', 0, 'bin_centers', [], ...
    'bin_offsets', [], 'bin_min', [], 'bin_max', [], 'has_model', false);

if numel(ref_rts) < 2 || numel(target_rts) < 2
    return;
end

% Fit global linear model with optional outlier removal
ref_rts = ref_rts(:);
target_rts = target_rts(:);
p = polyfit(ref_rts, target_rts, 1);

outlier_k = COptionUtils.get(options, 'outlier_k', 3);
outlier_method = COptionUtils.get(options, 'outlier_method', 'mad');
if outlier_k > 0
    rt_pred = p(1) .* ref_rts + p(2);
    residuals = target_rts - rt_pred;
    if strcmpi(outlier_method, 'mad')
        mad_res = median(abs(residuals - median(residuals)));
        if mad_res > 0
            keep = abs(residuals - median(residuals)) <= outlier_k * mad_res;
        else
            keep = true(size(residuals));
        end
    else
        sigma = std(residuals);
        if sigma > 0
            keep = abs(residuals) <= outlier_k * sigma;
        else
            keep = true(size(residuals));
        end
    end
    if sum(keep) >= 2
        p = polyfit(ref_rts(keep), target_rts(keep), 1);
    end
end

model.slope = p(1);
model.intercept = p(2);
model.has_model = true;

num_bins = COptionUtils.get(options, 'num_bins', 5);
min_per_bin = COptionUtils.get(options, 'min_per_bin', 5);
if num_bins <= 1
    return;
end

% Calculate residuals for fitting local offsets in RT bins
rt_pred = model.slope .* ref_rts + model.intercept;
residuals = target_rts - rt_pred;

edges = linspace(min(ref_rts), max(ref_rts), num_bins + 1);
bin_centers = zeros(1, num_bins);
bin_offsets = zeros(1, num_bins);
has_bin = false(1, num_bins);

for idx = 1:num_bins
    bin_centers(idx) = mean([edges(idx), edges(idx+1)]);
    in_bin = ref_rts >= edges(idx) & ref_rts <= edges(idx+1);
    if sum(in_bin) >= min_per_bin
        bin_offsets(idx) = median(residuals(in_bin));
        has_bin(idx) = true;
    end
end

if any(has_bin)
    model.bin_centers = bin_centers(has_bin);
    model.bin_offsets = bin_offsets(has_bin);
    model.bin_min = min(model.bin_centers);
    model.bin_max = max(model.bin_centers);
end
end

