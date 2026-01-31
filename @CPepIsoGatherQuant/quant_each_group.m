function [has_nonzero_imp, imp_idx_nonzero, area_imp_final, rt_bound, idx_selected, ratio_each_XIC_peak]...
    = quant_each_group(obj,raw_name,ratio_raw,rt_raw,...
    intensity_raw, low_mz_bound, high_mz_bound, selected_charge)
% Quantify each group
% input:
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
% output:
%   has_nonzero_imp (1 x 1 logical)
%       is there non zero area under XIC
%   imp_idx_nonzero (M x 1 double)
%       indices of non-zero area IMPs (subset of 1..K)
%   area_imp_final (M x 1 double) area
%       total quantification of each selected IMP, area under curve of XIC
%   rt_bound (M x 1 struct)
%       retention time bound for each selected IMP, fields: .start/.end (minutes)
%   idx_selected (M x 1 double)
%       selected XIC peak index per selected IMP
%   ratio_each_XIC_peak (M x P double)
%       area ratio of each IMP within each candidate peak

has_nonzero_imp = false;
imp_idx_nonzero = [];
area_imp_final = [];
rt_bound = [];
idx_selected = [];
ratio_each_XIC_peak = [];

% Preprocess inputs (Sort, Smooth, Denoise)
[rt_sorted, ratio_sorted, xic_rt, xic_intensity_smoothed, xic_intensity_raw, is_valid] = ...
    CQuantIMPGroupUtils.prepare_ms1_xic(...
        obj.m_cMs12DatasetIO, raw_name, rt_raw, intensity_raw, ratio_raw, ...
        obj.m_minMSMSnum, low_mz_bound, high_mz_bound, selected_charge);

if ~is_valid
    return;
end

% Extract the XIC peaks around the identified MSMS precursor
XIC_peaks = CChromatogramUtils.detect_xic_peaks(xic_rt, xic_intensity_smoothed, xic_intensity_raw, rt_sorted, obj.m_alpha);

if isempty(XIC_peaks)
    return;
end

% Calculate the ratio on each XIC points using kernel method
ratio_estimated = CChromatogramUtils.calculate_kernel_ratio(xic_rt, rt_sorted, ratio_sorted, XIC_peaks, true);

% Peak-wise filtering and normalization
ratio_estimated = CQuantIMPGroupUtils.filter_and_normalize_peak_ratios(...
    xic_rt, xic_intensity_smoothed, ratio_estimated, XIC_peaks, obj.m_resFilterThres);

% For each IMP, evaluate all candidate XIC peaks by computing:
%   - imp_max_props, max peak contribution ratio
%   - peak_fwhms: half maximum peak width (not used currently)
%   - area_imp_by_peak: area contribution in each peak
% [imp_max_props, peak_fwhms, area_imp_by_peak, rt_bound] = ...
[imp_max_props, ~, area_imp_by_peak, rt_bound] = ...
    CQuantIMPGroupUtils.compute_peak_features(xic_rt, xic_intensity_smoothed, ratio_estimated, XIC_peaks);

% Peak Selection (Per IMP)
idx_selected = CQuantIMPGroupUtils.select_best_peak_per_imp(imp_max_props, area_imp_by_peak);

% Global Refinement (Re-distribution based on Selection)
% ratio_estimated = CQuantIMPGroupUtils.refine_ratios_by_selection(ratio_estimated, XIC_peaks, idx_selected);

% Final Area Calculation (reuse cached peak areas)
area_imp_final = CQuantIMPGroupUtils.get_final_area_from_peak_areas(...
    area_imp_by_peak, idx_selected);

% Get the non-zero area under XIC, index and rt_bound
[imp_idx_nonzero, area_imp_final, rt_bound, idx_selected, area_imp_by_peak] = ...
    CQuantIMPGroupUtils.filter_nonzero_xic(area_imp_final, rt_bound, ...
    idx_selected, area_imp_by_peak);
ratio_each_XIC_peak = area_imp_by_peak./sum(area_imp_by_peak,1);
if ~isempty(imp_idx_nonzero)
    has_nonzero_imp = true;
end
end