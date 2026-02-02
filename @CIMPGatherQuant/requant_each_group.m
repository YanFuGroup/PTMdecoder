function [has_nonzero_imp, imp_idx_nonzero, area_imp_final, rt_bound, max_label, ratio_each_XIC_peak]...
    = requant_each_group(obj,raw_name,ratio_raw,rt_raw,...
    intensity_raw, low_mz_bound, high_mz_bound, selected_charge,...
    current_imp_rt_range)
% Re-quantify each group
% input:
%   raw_name (1 x 1 char/string)
%       the name of the raw (mgf) file
%   ratio_raw (N x K double)
%       ratio matrix of quantification in current group; rows aligned to rt_raw
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
%   current_imp_rt_range (K x 1 cell)
%       RT ranges for each IMP (each cell: [] or [start, end] in minutes)
% output:
%   has_nonzero_imp (1 x 1 logical)
%       is there non zero area under XIC
%   imp_idx_nonzero (M x 1 double)
%       indices of non-zero area IMPs (subset of 1..K)
%   area_imp_final (M x 1 double) area
%       total quantification of each selected IMP, area under curve of XIC
%   rt_bound (M x 1 struct)
%       retention time bound, fields: .start/.end (minutes)
%   max_label (M x 1 double)
%       max check label of the XIC peak for each IMP
%   ratio_each_XIC_peak (M x 1 double)
%       ratio of each IMP area to total XIC area in its peak

has_nonzero_imp = false;
rt_error_tol = 1; % RT tolerance in minutes

% Preprocess inputs (Sort, Smooth, Denoise) and get Smoothed XIC
[rt_sorted, ratio_sorted, xic_rt, xic_intensity_smoothed, ~, is_valid] = ...
    CQuantIMPGroupPreprocessUtils.prepare_ms1_xic(...
        obj.m_cMs12DatasetIO, raw_name, rt_raw, intensity_raw, ratio_raw, ...
        obj.m_minMSMSnum, low_mz_bound, high_mz_bound, selected_charge);

if ~is_valid
    imp_idx_nonzero = [];
    area_imp_final = [];
    rt_bound = [];
    max_label = [];
    ratio_each_XIC_peak = [];
    return;
end

% Extract the rt bound of XIC peak and convert to index bounds
[final_XIC_peak_for_IMP, max_label, is_skip_vec, peak_ranges] = ...
    CQuantIMPGroupPreprocessUtils.prepare_peak_ranges_from_imp_rt_range(...
        xic_rt, current_imp_rt_range, rt_error_tol);

% Calculate the ratio on each XIC points using kernel method
ratio_estimated = CChromatogramUtils.calculate_kernel_ratio(xic_rt, rt_sorted, ratio_sorted, peak_ranges, false);


% Requantification using revised RT
[area_imp_final, rt_bound, ratio_each_XIC_peak] = ...
    CQuantIMPGroupAreaUtils.compute_imp_peak_area_and_ratio(...
        xic_rt, xic_intensity_smoothed, ratio_estimated, ...
        peak_ranges, final_XIC_peak_for_IMP, is_skip_vec);

% Get the non-zero area under XIC, index and rt_bound
[imp_idx_nonzero, area_imp_final, rt_bound, max_label, ratio_each_XIC_peak] = ...
    CQuantIMPGroupAreaUtils.filter_nonzero_xic(area_imp_final, rt_bound, max_label, ratio_each_XIC_peak);
if ~isempty(imp_idx_nonzero)
    has_nonzero_imp = true;
end
end

