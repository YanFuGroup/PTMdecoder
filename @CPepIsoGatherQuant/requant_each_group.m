function [bhave_non_zeros, idxNonZero, auxic, rt_bound, max_label, ratio_each_XIC_peak]...
    = requant_each_group(obj,raw_name,current_ratioMatrix,current_rts,...
    current_inten, low_mz_bound, high_mz_bound, selected_charge,...
    current_iso_rt_range)
% Re-quantify each group
% input:
%   raw_name (1 x 1 char/string)
%       the name of the raw (mgf) file
%   current_ratioMatrix (N x K double)
%       ratio matrix of quantification in current group; rows aligned to current_rts
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
%   current_iso_rt_range (K x 1 cell)
%       RT ranges for each IMP (each cell: [] or [start, end] in minutes)
% output:
%   bhave_non_zeros (1 x 1 logical)
%       is there non zero area under XIC
%   idxNonZero (M x 1 double)
%       indices of non-zero area IMPs (subset of 1..K)
%   auxic (M x 1 double) area
%       total quantification of each selected IMP, area under curve of XIC
%   rt_bound (M x 1 struct)
%       retention time bound, fields: .start/.end (minutes)
%   max_label (M x 1 double)
%       max check label of the XIC peak for each IMP
%   ratio_each_XIC_peak (M x 1 double)
%       ratio of each IMP area to total XIC area in its peak

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
    CQuantIMPGroupUtils.prepare_peak_ranges_from_imp_rt_range(...
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

