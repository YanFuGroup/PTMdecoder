function [has_nonzero_imp, imp_idx_nonzero, area_imp_final, xic_peak_rt_bounds, idx_selected, ratio_each_XIC_peak] = ...
    quantGroup(cMs12DatasetIO, raw_name, ratio_raw, rt_raw, intensity_raw, low_mz_bound, high_mz_bound, selected_charge, minMSMSnum, alpha, resFilterThres)
% Quantify each group
% input:
%   cMs12DatasetIO (object)
%       dataset IO object
%   raw_name (1 x 1 char/string)
%       the name of the raw (mgf) file
%   ratio_raw (N x K double)
%       ratio matrix of quantification in current group; rows aligned to rt_raw
%       and columns are IMPs; each row sums to ~1 within a peak after normalization
%   rt_raw (N x 1 double) minutes
%       retention time in current group
%   intensity_raw (N x 1 double) intensity
%       intensity in current group
%   low_mz_bound (1 x 1 double) m/z
%       low precursor m/z bound
%   high_mz_bound (1 x 1 double) m/z
%       high precursor m/z bound
%   selected_charge (1 x 1 double/int)
%       current precursor charge
%   minMSMSnum (1 x 1 double/int)
%       minimum MSMS number
%   alpha (1 x 1 double)
%       threshold parameter
%   resFilterThres (1 x 1 double)
%       resolution filter threshold
% output:
%   has_nonzero_imp (1 x 1 logical)
%       is there non zero area under XIC
%   imp_idx_nonzero (M x 1 double)
%       indices of non-zero area IMPs (subset of 1..K)
%   area_imp_final (M x 1 double) area
%       total quantification of each selected IMP, area under curve of XIC
%   xic_peak_rt_bounds (struct array)
%       retention time bounds for peaks
%   idx_selected (M x 1 double)
%       indices of selected peaks
%   ratio_each_XIC_peak (M x P double)
%       ratio for each XIC peak

has_nonzero_imp = false;
imp_idx_nonzero = [];
area_imp_final = [];
xic_peak_rt_bounds = [];
idx_selected = [];
ratio_each_XIC_peak = [];

% Preprocess inputs (Sort, Smooth, Denoise)
[rt_sorted, ratio_sorted, xic_rt, xic_intensity_smoothed, xic_intensity_raw, is_valid] = ...
    CXICPreprocessUtils.prepare_ms1_xic(...
        cMs12DatasetIO, raw_name, rt_raw, intensity_raw, ratio_raw, ...
        minMSMSnum, low_mz_bound, high_mz_bound, selected_charge);

if ~is_valid
    return;
end

% Record rt_sorted length for offline analysis
CIMPQuantifier.rt_sorted_stats('record', numel(rt_sorted));

% Extract the XIC peaks around the identified MSMS precursor
xic_peak_idx_bounds = CXICSignalUtils.detect_xic_peaks(xic_rt, xic_intensity_smoothed, xic_intensity_raw, rt_sorted, alpha);

if isempty(xic_peak_idx_bounds)
    return;
end

% Calculate the ratio on each XIC point using kernel method
xic_ratio_estimated = CXICPeakUtils.calculate_kernel_ratio(xic_rt, rt_sorted, ratio_sorted, xic_peak_idx_bounds, true);

% Peak-wise filtering and normalization
xic_ratio_estimated = CXICPeakUtils.filter_and_normalize_peak_ratios(...
    xic_rt, xic_intensity_smoothed, xic_ratio_estimated, xic_peak_idx_bounds, resFilterThres);

% For each IMP, evaluate all candidate XIC peaks by computing:
%   - imp_max_props, max peak contribution ratio
%   - peak_fwhms: half maximum peak width (not used currently)
%   - area_imp_by_peak: area contribution in each peak
% [imp_max_props, peak_fwhms, area_imp_by_peak, xic_peak_rt_bounds] = ...
[imp_max_props, ~, area_imp_by_peak, xic_peak_rt_bounds] = ...
    CXICPeakUtils.compute_peak_features(xic_rt, xic_intensity_smoothed, xic_ratio_estimated, xic_peak_idx_bounds);

% Peak Selection (Per IMP)
idx_selected = CXICPeakUtils.select_best_peak_per_imp(imp_max_props, area_imp_by_peak);

% Global Refinement (Re-distribution based on Selection)
% xic_ratio_estimated = CXICPeakUtils.refine_ratios_by_selection(xic_ratio_estimated, xic_peak_idx_bounds, idx_selected);

% Final Area Calculation (reuse cached peak areas)
area_imp_final = CXICAreaUtils.get_final_area_from_peak_areas(...
    area_imp_by_peak, idx_selected);

% Get the non-zero area under XIC, index and xic_peak_rt_bounds
[imp_idx_nonzero, area_imp_final, xic_peak_rt_bounds, idx_selected, area_imp_by_peak] = ...
    CXICAreaUtils.filter_nonzero_xic(area_imp_final, xic_peak_rt_bounds, ...
    idx_selected, area_imp_by_peak);
sums = sum(area_imp_by_peak, 1);
ratio_each_XIC_peak = zeros(size(area_imp_by_peak));
valid_peaks = sums > 0;
ratio_each_XIC_peak(:, valid_peaks) = area_imp_by_peak(:, valid_peaks) ./ sums(valid_peaks);
if ~isempty(imp_idx_nonzero)
    has_nonzero_imp = true;
end
end
