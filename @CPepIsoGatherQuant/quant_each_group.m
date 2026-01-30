function [bhave_non_zeros, idxNonZero, auxic, rt_bound, idx_selected, ratio_each_XIC_peak]...
    = quant_each_group(obj,raw_name,current_ratioMatrix,current_rts,...
    current_inten, low_mz_bound, high_mz_bound, selected_charge)
% Quantify each group
% input:    
%   raw_name
%       the name of the raw (mgf) file
%   current_iso_name
%       names of IMPs in current group
%   current_ratioMatrix
%       ratio matrix of quantification in current group
%   current_rts
%       retention time in current group
%   current_inten
%       intensity in current group
%   low_mz_bound
%       low precursor m/z bound
%   high_mz_bound
%       high precursor m/z bound
%   selected_charge
%       current precursor charge
% output:   
%   bhave_non_zeros
%       is there non zero area under XIC
%   idxNonZero
%       the indices of non zero area under XIC
%   area
%       total quantification of each IMP in current group,
%       area under curve of XIC
%   rt_bound
%       the retention time bound, .start and .end
%   idx_selected
%       the indices of finally selected XIC peak

bhave_non_zeros = false;
idxNonZero = [];
auxic = [];
rt_bound = [];
idx_selected = [];
ratio_each_XIC_peak = [];
num_imp = size(current_ratioMatrix,2);

% Preprocess inputs (Sort, Smooth, Denoise)
[sort_rts, sort_inten, sort_ratioMatrix, is_valid] = ...
    CChromatogramUtils.preprocess_ms1_inputs(current_rts, current_inten, current_ratioMatrix, obj.m_minMSMSnum);

if ~is_valid
    bhave_non_zeros = false;
    idxNonZero = [];
    auxic = [];
    rt_bound = [];
    idx_selected = [];
    ratio_each_XIC_peak = [];
    return;
end

% Get Smoothed XIC
[rt_grid, smoothed_intensity, intensity] = ...
    CChromatogramUtils.get_smoothed_xic(obj.m_cMs12DatasetIO, raw_name, low_mz_bound, high_mz_bound, selected_charge);

% Extract the XIC peaks around the identified MSMS precursor
XIC_peaks = CChromatogramUtils.detect_xic_peaks(rt_grid, smoothed_intensity, intensity, sort_rts, obj.m_alpha);

if isempty(XIC_peaks)
    return;
end

% Calculate the ratio on each XIC points using kernel method
esti_ratio = CChromatogramUtils.calculate_kernel_ratio(rt_grid, sort_rts, sort_ratioMatrix, XIC_peaks, true);

% Stage 1: peak-wise filtering and normalization
esti_ratio = CQuantIMPGroupUtils.filter_and_normalize_peak_ratios(...
    rt_grid, smoothed_intensity, esti_ratio, XIC_peaks, obj.m_resFilterThres);

% ========================================================================
% Stage 2: Feature Calculation (Vectorized over Peaks)
% ========================================================================

% For each IMP, evaluate all candidate XIC peaks by computing:
%   - imp_max_props, max peak contribution ratio
%   - peak_fwhms: half maximum peak width
%   - ratio_each_XIC_peak: area contribution in each peak
auxic = zeros(num_imp,1); 
idx_selected = zeros(num_imp,1);
[imp_max_props, peak_fwhms, ratio_each_XIC_peak, rt_bound] = ...
    CQuantIMPGroupUtils.compute_peak_features(rt_grid, smoothed_intensity, esti_ratio, XIC_peaks);

% ========================================================================
% Stage 3: Peak Selection (Per IMP)
% ========================================================================
idx_selected = CQuantIMPGroupUtils.select_best_peak_per_imp(imp_max_props, ratio_each_XIC_peak);

% ========================================================================
% Stage 4: Global Refinement (Re-distribution based on Selection)
% ========================================================================
esti_ratio = CQuantIMPGroupUtils.refine_ratios_by_selection(esti_ratio, XIC_peaks, idx_selected);

% ========================================================================
% Stage 5: Final Area Calculation
% ========================================================================
auxic = CQuantIMPGroupUtils.compute_final_area(...
    rt_grid, smoothed_intensity, esti_ratio, XIC_peaks, idx_selected);

% Get the non-zero area under XIC, index and rt_bound
idxNonZero = find(auxic(:,1)~=0);
auxic = auxic(idxNonZero,:);
rt_bound = rt_bound(idxNonZero,:);
idx_selected = idx_selected(idxNonZero,:);
% TODO: use auxic to calculate the ratio_each_XIC_peak instead.
ratio_each_XIC_peak = ratio_each_XIC_peak(idxNonZero,:);
ratio_each_XIC_peak = ratio_each_XIC_peak./sum(ratio_each_XIC_peak,1);
if ~isempty(idxNonZero)
    bhave_non_zeros = true;
end
end