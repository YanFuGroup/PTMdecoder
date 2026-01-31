function [area_imp_final, rt_bound, ratio_each_XIC_peak] = compute_imp_peak_area_and_ratio(...
    xic_rt, xic_intensity_smoothed, ratio_estimated, peak_ranges, final_XIC_peak_for_IMP, is_skip_vec)
% Compute area and ratio for each IMP based on peak ranges.
% input:
%   xic_rt (N x 1 double) minutes
%       retention time grid
%   xic_intensity_smoothed (N x 1 double) intensity
%       total smoothed XIC intensity
%   ratio_estimated (N x K double)
%       estimated ratio of each IMP across RT grid
%   peak_ranges (K x 1 struct)
%       index bounds for each IMP peak; fields: left_bound/right_bound (indices into xic_rt)
%   final_XIC_peak_for_IMP (K x 1 struct)
%       RT bounds for each IMP peak; fields: left_bound/right_bound (minutes)
%   is_skip_vec (K x 1 logical)
%       vector indicating IMPs to skip
% output:
%   area_imp_final (K x 1 double) area
%       area under each IMP XIC
%   rt_bound (K x 1 struct)
%       RT bounds for each IMP; fields: .start/.end (minutes)
%   ratio_each_XIC_peak (K x 1 double)
%       ratio of each IMP area to total XIC area in its peak

intensityMatrix = ratio_estimated.*xic_intensity_smoothed;
num_imp = size(intensityMatrix, 2);
area_imp_final = zeros(num_imp,1);
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
    area_imp_final(idx_imp,1) = CChromatogramUtils.calculate_area(...
        xic_rt, intensityMatrix(:,idx_imp), idx_start, idx_end);

    % Calculate total area for ratio
    % NOTE: Potential optimization if peak_ranges have many repeats:
    % cache total_temp for unique (idx_start, idx_end) pairs to reduce
    % repeated calculate_area calls on xic_intensity_smoothed.
    total_temp = CChromatogramUtils.calculate_area(...
        xic_rt, xic_intensity_smoothed, idx_start, idx_end);

    ratio_each_XIC_peak(idx_imp,1) = area_imp_final(idx_imp,1) / total_temp;
end
end
