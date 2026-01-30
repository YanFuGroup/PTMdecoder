function [bhave_non_zeros, idxNonZero, auxic, rt_bound, idx_selected, ratio_each_XIC_peak]...
    = quant_each_group(obj,raw_name,current_ratioMatrix,current_rts,...
    current_inten, low_mz_bound, high_mz_bound, selected_charge)
% Quantify each group
% input:
%   raw_name (1 x 1 char/string)
%       the name of the raw (mgf) file
%   current_ratioMatrix (N x K double)
%       ratio matrix of quantification in current group; rows aligned to current_rts
%       and columns are IMPs; each row sums to ~1 within a peak after normalization
%   current_rts (N x 1 double) minutes
%       retention time in current group
%   current_inten (N x 1 double) intensity
%       intensity in current group
%   low_mz_bound (1 x 1 double) m/z
%       low precursor m/z bound
%   high_mz_bound (1 x 1 double) m/z
%       high precursor m/z bound
%   selected_charge (1 x 1 double/int)
%       current precursor charge
% output:
%   bhave_non_zeros (1 x 1 logical)
%       is there non zero area under XIC
%   idxNonZero (M x 1 double)
%       indices of non-zero area IMPs (subset of 1..K)
%   auxic (M x 1 double) area
%       total quantification of each selected IMP, area under curve of XIC
%   rt_bound (M x 1 struct)
%       retention time bound for each selected IMP, fields: .start/.end (minutes)
%   idx_selected (M x 1 double)
%       selected XIC peak index per selected IMP
%   ratio_each_XIC_peak (M x P double)
%       area ratio of each IMP within each candidate peak

bhave_non_zeros = false;
idxNonZero = [];
auxic = [];
rt_bound = [];
idx_selected = [];
ratio_each_XIC_peak = [];

% Preprocess inputs (Sort, Smooth, Denoise)
[sort_rts, sort_ratioMatrix, rt_grid, smoothed_intensity, intensity, is_valid] = ...
    CQuantIMPGroupUtils.prepare_ms1_xic(...
        obj.m_cMs12DatasetIO, raw_name, current_rts, current_inten, current_ratioMatrix, ...
        obj.m_minMSMSnum, low_mz_bound, high_mz_bound, selected_charge);

if ~is_valid
    return;
end

% Extract the XIC peaks around the identified MSMS precursor
XIC_peaks = CChromatogramUtils.detect_xic_peaks(rt_grid, smoothed_intensity, intensity, sort_rts, obj.m_alpha);

if isempty(XIC_peaks)
    return;
end

% Calculate the ratio on each XIC points using kernel method
esti_ratio = CChromatogramUtils.calculate_kernel_ratio(rt_grid, sort_rts, sort_ratioMatrix, XIC_peaks, true);

% Peak-wise filtering and normalization
esti_ratio = CQuantIMPGroupUtils.filter_and_normalize_peak_ratios(...
    rt_grid, smoothed_intensity, esti_ratio, XIC_peaks, obj.m_resFilterThres);

% For each IMP, evaluate all candidate XIC peaks by computing:
%   - imp_max_props, max peak contribution ratio
%   - peak_fwhms: half maximum peak width (not used currently)
%   - area_each_XIC_peak: area contribution in each peak
% [imp_max_props, peak_fwhms, area_each_XIC_peak, rt_bound] = ...
[imp_max_props, ~, area_each_XIC_peak, rt_bound] = ...
    CQuantIMPGroupUtils.compute_peak_features(rt_grid, smoothed_intensity, esti_ratio, XIC_peaks);

% Peak Selection (Per IMP)
idx_selected = CQuantIMPGroupUtils.select_best_peak_per_imp(imp_max_props, area_each_XIC_peak);

% Global Refinement (Re-distribution based on Selection)
% esti_ratio = CQuantIMPGroupUtils.refine_ratios_by_selection(esti_ratio, XIC_peaks, idx_selected);

% Final Area Calculation (reuse cached peak areas)
auxic = CQuantIMPGroupUtils.get_final_area_from_peak_areas(...
    area_each_XIC_peak, idx_selected);

% Get the non-zero area under XIC, index and rt_bound
[idxNonZero, auxic, rt_bound, idx_selected, area_each_XIC_peak] = ...
    CQuantIMPGroupUtils.filter_nonzero_xic(auxic, rt_bound, idx_selected, area_each_XIC_peak);
% TODO: use auxic to calculate the area_each_XIC_peak instead.
ratio_each_XIC_peak = area_each_XIC_peak./sum(area_each_XIC_peak,1);
if ~isempty(idxNonZero)
    bhave_non_zeros = true;
end
end