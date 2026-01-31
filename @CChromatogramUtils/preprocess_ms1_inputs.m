function [rt_sorted, ratio_sorted, is_valid] = preprocess_ms1_inputs(rt_raw, intensity_raw, ratio_raw, minMSMSnum)
    % Preprocess MS1 inputs: Sort by retention time, Smooth intensity, and Denoise
    % Inputs:
    %   rt_raw (N x 1 double) minutes
    %       Array of retention times
    %   intensity_raw (N x 1 double) intensity
    %       Array of intensities
    %   ratio_raw (N x K double)
    %       Quantification matrix
    %   minMSMSnum (1 x 1 double/int)
    %       Minimum number of MSMS recurrence required
    % Output:
    %   rt_sorted (M x 1 double) minutes
    %       Sorted retention times after filtering
    %   ratio_sorted (M x K double)
    %       Sorted ratio matrix after filtering
    %   is_valid (1 x 1 logical)
    %       Whether inputs pass minimum MSMS count
    
    is_valid = true;
    
    % Sort MS1 signal (pair of retention time and intensity) by time
    rt_sorted = [(1:length(rt_raw))',rt_raw];
    rt_sorted = sortrows(rt_sorted,2); % Sort in ascending order by retention time column
    sort_idx = rt_sorted(:,1);
    rt_sorted = rt_sorted(:,2);
    intensity_sorted = intensity_raw(sort_idx);
    ratio_sorted = ratio_raw(sort_idx,:); % Rearrange the matrix in chronological order
    
    % Sort and denoise using a relative abundance threshold method
    maxInten = max(intensity_sorted);
    if isempty(maxInten)
        maxInten = 0;
    end
    
    tmp = intensity_sorted < 0.05 * maxInten; % Find results where intensity is less than 0.05 of the maximum abundance and discard them
    % intensity_sorted(tmp) = []; 
    rt_sorted(tmp) = [];
    ratio_sorted(tmp,:) = [];

    % Check if there are enough rows left
    % Replaces obj.hasMinRows(ratio_sorted, minMSMSnum)
    if size(ratio_sorted, 1) < minMSMSnum
        % If the ratio matrix has less than min rows, skip this group
        is_valid = false;
        ratio_sorted = [];
    end
    % TODO: The filtering is not only here, but also the XIC integration area is determined later.
    %  Then, within the integration area, it is checked how many fall within this area.
    %  If there are not enough, they will be removed.
end
