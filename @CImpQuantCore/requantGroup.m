function [has_nonzero_imp, imp_idx_nonzero, area_imp_final, xic_peak_rt_bounds, max_label, ratio_each_XIC_peak] = ...
    requantGroup(cMs12DatasetIO, raw_name, ratio_raw, rt_raw, intensity_raw, low_mz_bound, high_mz_bound, selected_charge, current_imp_rt_range, minMSMSnum)
% Re-quantify each group
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
%   current_imp_rt_range (array)
%       current IMP retention time range
%   minMSMSnum (1 x 1 double/int)
%       minimum MSMS number
% output:
%   has_nonzero_imp (1 x 1 logical)
%       is there non zero area under XIC
%   imp_idx_nonzero (M x 1 double)
%       indices of non-zero area IMPs (subset of 1..K)
%   area_imp_final (M x 1 double) area
%       total quantification of each selected IMP, area under curve of XIC
%   xic_peak_rt_bounds (struct array)
%       retention time bounds for peaks
%   max_label (array)
%       max labels
%   ratio_each_XIC_peak (M x P double)
%       ratio for each XIC peak

has_nonzero_imp = false;
imp_idx_nonzero = [];
area_imp_final = [];
xic_peak_rt_bounds = [];
max_label = [];
ratio_each_XIC_peak = [];

% Preprocess inputs (Sort, Smooth, Denoise)
[rt_sorted, ratio_sorted, xic_rt, xic_intensity_smoothed, ~, is_valid] = ...
    CXICPreprocessUtils.prepare_ms1_xic(...
        cMs12DatasetIO, raw_name, rt_raw, intensity_raw, ratio_raw, ...
        minMSMSnum, low_mz_bound, high_mz_bound, selected_charge);

if ~is_valid
    return;
end

% Filter peaks by manual RT range
[max_label, is_skip_vec, xic_peak_idx_bounds] = CXICPreprocessUtils.prepare_peak_ranges_from_imp_rt_range(...
    xic_rt, current_imp_rt_range, 1);

if isempty(xic_peak_idx_bounds)
    return;
end

% Calculate the ratio on each XIC point using kernel method
xic_ratio_estimated = CXICPeakUtils.calculate_kernel_ratio(xic_rt, rt_sorted, ratio_sorted, xic_peak_idx_bounds, false);

% Requantification using revised RT
[area_imp_final, xic_peak_rt_bounds, ratio_each_XIC_peak] = ...
    CXICAreaUtils.compute_imp_peak_area_and_ratio(...
        xic_rt, xic_intensity_smoothed, xic_ratio_estimated, ...
        xic_peak_idx_bounds, is_skip_vec);

% Get the non-zero area under XIC, index and xic_peak_rt_bounds
[imp_idx_nonzero, area_imp_final, xic_peak_rt_bounds, max_label, ratio_each_XIC_peak] = ...
    CXICAreaUtils.filter_nonzero_xic(area_imp_final, xic_peak_rt_bounds, max_label, ratio_each_XIC_peak);
if ~isempty(imp_idx_nonzero)
    has_nonzero_imp = true;
end
end
