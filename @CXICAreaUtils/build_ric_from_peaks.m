function ric = build_ric_from_peaks(xic_rt, xic_intensity_smoothed, xic_ratio_estimated, xic_peak_idx_bounds, is_skip_vec)
% Build RIC cell array for each IMP using closed peak data.
% input:
%   xic_rt (N x 1 double) minutes
%       retention time grid
%   xic_intensity_smoothed (N x 1 double) intensity
%       total smoothed XIC intensity
%   xic_ratio_estimated (N x K double)
%       estimated ratio of each IMP across RT grid
%   xic_peak_idx_bounds (K x 1 struct)
%       index bounds for each IMP peak; fields: idx_start/idx_end (indices into xic_rt)
%   is_skip_vec (K x 1 logical)
%       vector indicating IMPs to skip
% output:
%   ric (K x 2 cell)
%       cell array with rt and intensity per IMP; ric{i,1}=rt (minutes), ric{i,2}=intensity

intensityMatrix = xic_ratio_estimated.*xic_intensity_smoothed;
num_imp = size(intensityMatrix, 2);
ric = cell(num_imp, 2);
for idx_imp = 1:num_imp
    % Check if need to consider this IMP
    if is_skip_vec(idx_imp)
        continue;
    end

    % Retrieve closed peak data for plotting/integration validation
    [rec_rt, rec_inten] = CXICSignalUtils.get_closed_peak_data(...
        xic_rt, intensityMatrix(:,idx_imp), ...
        xic_peak_idx_bounds(idx_imp).idx_start, xic_peak_idx_bounds(idx_imp).idx_end);

    ric{idx_imp,1} = rec_rt;
    ric{idx_imp,2} = rec_inten;
end
end
