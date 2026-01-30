function [bhave_non_zeros, idxNonZero, auxic, rt_bound, max_label, ratio_each_XIC_peak]...
    = requant_each_group(obj,raw_name,current_ratioMatrix,current_rts,...
    current_inten, low_mz_bound, high_mz_bound, selected_charge,...
    current_iso_rt_range)
% Re-quantify each group
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
%   current_iso_rt_range
%       retention times of current IMPs
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
%   max_label
%       the max check label of the XIC peaks for each IMP

bhave_non_zeros = false;
rt_error_tol = 1; % RT tolerance in minutes

% Preprocess inputs (Sort, Smooth, Denoise) and get Smoothed XIC
[sort_rts, sort_ratioMatrix, rt_grid, smoothed_intensity, ~, is_valid] = ...
    CQuantIMPGroupUtils.prepare_ms1_xic(...
        obj.m_cMs12DatasetIO, raw_name, current_rts, current_inten, current_ratioMatrix, ...
        obj.m_minMSMSnum, low_mz_bound, high_mz_bound, selected_charge);

if ~is_valid
    idxNonZero = [];
    auxic = [];
    rt_bound = [];
    max_label = [];
    ratio_each_XIC_peak = [];
    return;
end

% Extract the rt bound of XIC peak and convert to index bounds
[final_XIC_peak_for_IMP, max_label, is_skip_vec, peak_ranges] = ...
    CQuantIMPGroupUtils.prepare_peak_ranges_from_iso_rt_range(...
        rt_grid, current_iso_rt_range, rt_error_tol);

% Calculate the ratio on each XIC points using kernel method
esti_ratio = CChromatogramUtils.calculate_kernel_ratio(rt_grid, sort_rts, sort_ratioMatrix, peak_ranges, false);


% Requantification using revised RT
[auxic, rt_bound, ratio_each_XIC_peak] = ...
    CQuantIMPGroupUtils.compute_imp_peak_area_and_ratio(...
        rt_grid, smoothed_intensity, esti_ratio, ...
        peak_ranges, final_XIC_peak_for_IMP, is_skip_vec);

% Get the non-zero area under XIC, index and rt_bound
[idxNonZero, auxic, rt_bound, max_label, ratio_each_XIC_peak] = ...
    CQuantIMPGroupUtils.filter_nonzero_xic(auxic, rt_bound, max_label, ratio_each_XIC_peak);
if ~isempty(idxNonZero)
    bhave_non_zeros = true;
end
end

