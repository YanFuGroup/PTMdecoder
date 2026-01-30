function [auxic, rt_bound, ratio_each_XIC_peak] = compute_imp_peak_area_and_ratio(...
    rt_grid, smoothed_intensity, esti_ratio, peak_ranges, final_XIC_peak_for_IMP, is_skip_vec)
% Compute area and ratio for each IMP based on peak ranges.
% input:
%   rt_grid (N x 1 double) minutes
%       retention time grid
%   smoothed_intensity (N x 1 double) intensity
%       total smoothed XIC intensity
%   esti_ratio (N x K double)
%       estimated ratio of each IMP across RT grid
%   peak_ranges (K x 1 struct)
%       index bounds for each IMP peak; fields: left_bound/right_bound (indices into rt_grid)
%   final_XIC_peak_for_IMP (K x 1 struct)
%       RT bounds for each IMP peak; fields: left_bound/right_bound (minutes)
%   is_skip_vec (K x 1 logical)
%       vector indicating IMPs to skip
% output:
%   auxic (K x 1 double) area
%       area under each IMP XIC
%   rt_bound (K x 1 struct)
%       RT bounds for each IMP; fields: .start/.end (minutes)
%   ratio_each_XIC_peak (K x 1 double)
%       ratio of each IMP area to total XIC area in its peak

intensityMatrix = esti_ratio.*smoothed_intensity;
num_imp = size(intensityMatrix, 2);
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
end
