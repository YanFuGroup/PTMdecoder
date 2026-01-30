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
num_imp = size(current_ratioMatrix,2);
rt_error_tol = 1; % RT tolerance in minutes
% A vector showing is needed to skip this IMP.
%   Cannot delete because other filter can also change this vector
is_skip_vec = cellfun(@isempty,current_iso_rt_range);

% Preprocess inputs (Sort, Smooth, Denoise)
[sort_rts, sort_ratioMatrix, is_valid] = ...
    CChromatogramUtils.preprocess_ms1_inputs(current_rts, current_inten, current_ratioMatrix, obj.m_minMSMSnum);

if ~is_valid
    idxNonZero = [];
    auxic = [];
    rt_bound = [];
    max_label = [];
    ratio_each_XIC_peak = [];
    return;
end

% Get Smoothed XIC
[rt_grid, smoothed_intensity, ~] = ...
    CChromatogramUtils.get_smoothed_xic(obj.m_cMs12DatasetIO, raw_name, low_mz_bound, high_mz_bound, selected_charge);

% Extract the rt bound of XIC peak
[final_XIC_peak_for_IMP, max_label, is_skip_vec] = ...
    CChromatogramUtils.parse_imp_rt_ranges(current_iso_rt_range, is_skip_vec);

% Convert RT bounds to index bounds for CChromatogramUtils
peak_ranges = CChromatogramUtils.map_rt_to_indices(rt_grid, final_XIC_peak_for_IMP, is_skip_vec, rt_error_tol);

% Calculate the ratio on each XIC points using kernel method
esti_ratio = CChromatogramUtils.calculate_kernel_ratio(rt_grid, sort_rts, sort_ratioMatrix, peak_ranges, false);


% Requantification using revised RT
intensityMatrix = esti_ratio.*smoothed_intensity;
auxic = zeros(num_imp,1);
rt_bound = repmat(struct('start',0,'end',0), num_imp, 1);
ratio_each_XIC_peak = zeros(num_imp,1);

for idx_imp = 1:num_imp
    % Check if need to consider this IMP
    if is_skip_vec(idx_imp)
        continue;
    end
    
    % Use pre-calculated indices from map_rt_to_indices
    idx_start = peak_ranges(idx_imp).left_bound;
    idx_end = peak_ranges(idx_imp).right_bound;

    % Assign rt_bound from input peaks
    rt_bound(idx_imp).start = final_XIC_peak_for_IMP(idx_imp).left_bound;
    rt_bound(idx_imp).end = final_XIC_peak_for_IMP(idx_imp).right_bound;

    % Calculate area using the closed peak logic
    auxic(idx_imp,1) = CChromatogramUtils.calculate_area(...
        rt_grid, intensityMatrix(:,idx_imp), idx_start, idx_end);

    % Calculate total area for ratio
    % NOTE: Potential optimization if peak_ranges have many repeats:
    % cache total_temp for unique (idx_start, idx_end) pairs to reduce
    % repeated calculate_area calls on smoothed_intensity.
    total_temp = CChromatogramUtils.calculate_area(...
        rt_grid, smoothed_intensity, idx_start, idx_end);

    ratio_each_XIC_peak(idx_imp,1) = auxic(idx_imp,1) / total_temp;
end

% Get the non-zero area under XIC, index and rt_bound
[idxNonZero, auxic, rt_bound, max_label, ratio_each_XIC_peak] = ...
    CQuantIMPGroupUtils.filter_nonzero_xic(auxic, rt_bound, max_label, ratio_each_XIC_peak);
if ~isempty(idxNonZero)
    bhave_non_zeros = true;
end
end

