function [rt_peak_a, rt_peak_b, stats] = alignImpPair(obj, ms12DatasetIO, raw_name_a, imp_info_a, ...
    raw_name_b, imp_info_b, model, options, minMSMSnum, alpha, resFilterThres)
% Jointly select a peak pair for one IMP across two runs.
% Input:
%   obj (CXICAligner)
%       Aligner instance
%   ms12DatasetIO (CMS12DatasetIO)
%       MS1/MS2 dataset IO
%   raw_name_a (1 x 1 char/string)
%       Raw name of run A
%   imp_info_a (struct)
%       IMP info for run A (impName, impMass, ratio, rts, intensity, charge, lowMzBound, highMzBound)
%   raw_name_b (1 x 1 char/string)
%       Raw name of run B
%   imp_info_b (struct)
%       IMP info for run B (impName, impMass, ratio, rts, intensity, charge, lowMzBound, highMzBound)
%   model (struct)
%       Alignment transform model (A -> B)
%   options (struct, optional)
%       - rt_sigma (1 x 1 double)
%           RT Gaussian sigma (minutes) for peak selection. Default: 0.5.
%       - max_rt_residual (1 x 1 double)
%           Max allowed RT residual (minutes) for peak pairing.
%       - dead_time_min (1 x 1 double)
%           Min allowed RT (minutes) for peak start; earlier peaks are ignored.
%   minMSMSnum (1 x 1 double/int, optional)
%       Minimum MSMS count for XIC preprocessing
%   alpha (1 x 1 double, optional)
%       Peak detection threshold factor
%   resFilterThres (1 x 1 double, optional)
%       Ratio filter threshold
% Output:
%   rt_peak_a (struct or [])
%       Selected peak info for run A
%   rt_peak_b (struct or [])
%       Selected peak info for run B
%   stats (struct)
%       Counters for aligned IMPs

if nargin < 8
    options = obj.m_options;
end
if nargin < 9 || isempty(minMSMSnum)
    minMSMSnum = 1;
end
if nargin < 10 || isempty(alpha)
    alpha = 0.01;
end
if nargin < 11 || isempty(resFilterThres)
    resFilterThres = 0.01;
end

[data_a, ok_a] = compute_imp_peaks(ms12DatasetIO, raw_name_a, imp_info_a, minMSMSnum, alpha, resFilterThres);
[data_b, ok_b] = compute_imp_peaks(ms12DatasetIO, raw_name_b, imp_info_b, minMSMSnum, alpha, resFilterThres);

rt_peak_a = [];
rt_peak_b = [];
stats = struct('num_aligned', 0, 'num_total', 1);

if ~ok_a || ~ok_b
    return;
end

sigma = CStructOptionUtils.get(options, 'rt_sigma', 0.5);
max_rt_residual = CStructOptionUtils.get(options, 'max_rt_residual', []);
dead_time_min = CStructOptionUtils.get(options, 'dead_time_min', []);

scores_a = data_a.area_imp_by_peak .* data_a.imp_max_props;
scores_b = data_b.area_imp_by_peak .* data_b.imp_max_props;
if all(scores_a == 0) || all(scores_b == 0)
    return;
end

best_score = -Inf;
best_a_idx = 0;
best_b_idx = 0;
for idx_a = 1:numel(scores_a)
    if scores_a(idx_a) <= 0
        continue;
    end
    a_peak = data_a.xic_peak_rt_bounds(idx_a);
    if ~isempty(dead_time_min) && a_peak.rt_start < dead_time_min
        continue;
    end
    a_center = mean([a_peak.rt_start, a_peak.rt_end]);
    if model.has_model
        pred_b = CXICAlignTransform.apply(model, a_center);
    else
        warn_no_model_once(raw_name_a, raw_name_b);
        pred_b = a_center;
    end

    for idx_b = 1:numel(scores_b)
        if scores_b(idx_b) <= 0
            continue;
        end
        b_peak = data_b.xic_peak_rt_bounds(idx_b);
        if ~isempty(dead_time_min) && b_peak.rt_start < dead_time_min
            continue;
        end
        b_center = mean([b_peak.rt_start, b_peak.rt_end]);
        residual = abs(b_center - pred_b);
        if isempty(max_rt_residual) && isfield(model, 'rt_residual_threshold')
            max_rt_residual = model.rt_residual_threshold;
        end
        if ~isempty(max_rt_residual) && residual > max_rt_residual
            continue;
        end
        weight = exp(-residual / sigma);
        score = scores_a(idx_a) * scores_b(idx_b) * weight;
        if score > best_score
            best_score = score;
            best_a_idx = idx_a;
            best_b_idx = idx_b;
        end
    end
end

if best_a_idx > 0 && best_b_idx > 0
    a_peak = data_a.xic_peak_rt_bounds(best_a_idx);
    b_peak = data_b.xic_peak_rt_bounds(best_b_idx);
    rt_peak_a = struct( ...
        'rt_start', a_peak.rt_start, ...
        'rt_end', a_peak.rt_end, ...
        'ratio', data_a.ratio_each_XIC_peak(best_a_idx), ...
        'check_label', 2);
    rt_peak_b = struct( ...
        'rt_start', b_peak.rt_start, ...
        'rt_end', b_peak.rt_end, ...
        'ratio', data_b.ratio_each_XIC_peak(best_b_idx), ...
        'check_label', 2);
    stats.num_aligned = 1;
end
end


function warn_no_model_once(raw_name_a, raw_name_b)
% Warn once if no alignment model is available for the pair.
% Input:
%   raw_name_a (1 x 1 char/string)
%       Raw name of run A
%   raw_name_b (1 x 1 char/string)
%       Raw name of run B

persistent warned_pairs
if isempty(warned_pairs)
    warned_pairs = containers.Map('KeyType', 'char', 'ValueType', 'logical');
end
key = sprintf('%s->%s', raw_name_a, raw_name_b);
if ~isKey(warned_pairs, key)
    warning('CXICAligner:NoAlignmentModel', ...
        'No alignment model for pair %s. Using raw RT centers without alignment.', key);
    warned_pairs(key) = true;
end
end


function [data, is_ok] = compute_imp_peaks(ms12DatasetIO, raw_name, imp_info, minMSMSnum, alpha, resFilterThres)
% Extract peak features for one IMP.
% Input:
%   ms12DatasetIO (CMS12DatasetIO)
%       MS1/MS2 dataset IO
%   raw_name (1 x 1 char/string)
%       Raw name for XIC extraction
%   imp_info (struct)
%       IMP info (impName, impMass, ratio, rts, intensity, charge, lowMzBound, highMzBound)
%   minMSMSnum (1 x 1 double/int)
%       Minimum MSMS count for XIC preprocessing
%   alpha (1 x 1 double)
%       Peak detection threshold factor
%   resFilterThres (1 x 1 double)
%       Ratio filter threshold
% Output:
%   data (struct)
%       Peak feature container for one IMP
%   is_ok (logical)
%       True if peaks are detected

data = struct('imp_max_props', [], 'area_imp_by_peak', [], ...
    'xic_peak_rt_bounds', [], 'ratio_each_XIC_peak', []);

[rt_sorted, ratio_sorted, xic_rt, xic_intensity_smoothed, xic_intensity_raw, is_valid] = ...
    CXICPreprocessUtils.prepare_ms1_xic(...
        ms12DatasetIO, raw_name, imp_info.rts, ...
        imp_info.intensity, imp_info.ratio, ...
        minMSMSnum, imp_info.lowMzBound, imp_info.highMzBound, imp_info.charge);
if ~is_valid
    is_ok = false;
    return;
end

xic_peak_idx_bounds = CXICSignalUtils.detect_xic_peaks(...
    xic_rt, xic_intensity_smoothed, xic_intensity_raw, rt_sorted, alpha);
if isempty(xic_peak_idx_bounds)
    is_ok = false;
    return;
end

xic_ratio_estimated = CXICPeakUtils.calculate_kernel_ratio(...
    xic_rt, rt_sorted, ratio_sorted, xic_peak_idx_bounds, true);
xic_ratio_estimated = CXICPeakUtils.filter_and_normalize_peak_ratios(...
    xic_rt, xic_intensity_smoothed, xic_ratio_estimated, xic_peak_idx_bounds, resFilterThres);

[imp_max_props, ~, area_imp_by_peak, xic_peak_rt_bounds] = ...
    CXICPeakUtils.compute_peak_features(...
        xic_rt, xic_intensity_smoothed, xic_ratio_estimated, xic_peak_idx_bounds);

if isempty(imp_max_props)
    is_ok = false;
    return;
end

sums = sum(area_imp_by_peak, 1);
ratio_each_XIC_peak = zeros(size(area_imp_by_peak));
valid_peaks = sums > 0;
ratio_each_XIC_peak(:, valid_peaks) = area_imp_by_peak(:, valid_peaks) ./ sums(valid_peaks);

is_ok = true;
data.imp_max_props = imp_max_props;
data.area_imp_by_peak = area_imp_by_peak;
data.xic_peak_rt_bounds = xic_peak_rt_bounds;
data.ratio_each_XIC_peak = ratio_each_XIC_peak;
end
