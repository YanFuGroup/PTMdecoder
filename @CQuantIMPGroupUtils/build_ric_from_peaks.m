function ric = build_ric_from_peaks(rt_grid, smoothed_intensity, esti_ratio, peak_ranges, is_skip_vec)
% Build RIC cell array for each IMP using closed peak data.
% input:
%   rt_grid (N x 1 double) minutes
%       retention time grid
%   smoothed_intensity (N x 1 double) intensity
%       total smoothed XIC intensity
%   esti_ratio (N x K double)
%       estimated ratio of each IMP across RT grid
%   peak_ranges (K x 1 struct)
%       index bounds for each IMP peak; fields: left_bound/right_bound (indices into rt_grid)
%   is_skip_vec (K x 1 logical)
%       vector indicating IMPs to skip
% output:
%   ric (K x 2 cell)
%       cell array with rt and intensity per IMP; ric{i,1}=rt (minutes), ric{i,2}=intensity

intensityMatrix = esti_ratio.*smoothed_intensity;
num_imp = size(intensityMatrix, 2);
ric = cell(num_imp, 2);
for idx_imp = 1:num_imp
    % Check if need to consider this IMP
    if is_skip_vec(idx_imp)
        continue;
    end

    % Retrieve closed peak data for plotting/integration validation
    [rec_rt, rec_inten] = CChromatogramUtils.get_closed_peak_data(...
        rt_grid, intensityMatrix(:,idx_imp), ...
        peak_ranges(idx_imp).left_bound, peak_ranges(idx_imp).right_bound);

    ric{idx_imp,1} = rec_rt;
    ric{idx_imp,2} = rec_inten;
end
end
