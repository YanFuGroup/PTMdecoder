function [sort_rts, sort_ratioMatrix, is_valid] = preprocess_ms1_inputs(current_rts, current_inten, current_ratioMatrix, minMSMSnum)
    % Preprocess MS1 inputs: Sort by retention time, Smooth intensity, and Denoise
    % Inputs:
    %   current_rts: Array of retention times
    %   current_inten: Array of intensities
    %   current_ratioMatrix: quantification matrix
    %   minMSMSnum: Minimum number of MSMS recurrence required
    
    is_valid = true;
    
    % Sort MS1 signal (pair of retention time and intensity) by time
    sort_rts = [(1:length(current_rts))',current_rts];
    sort_rts = sortrows(sort_rts,2); % Sort in ascending order by retention time column
    sort_idx = sort_rts(:,1);
    sort_rts = sort_rts(:,2);
    sort_inten = current_inten(sort_idx);
    sort_ratioMatrix = current_ratioMatrix(sort_idx,:); % Rearrange the matrix in chronological order
    
    % Sort and denoise using a relative abundance threshold method
    maxInten = max(sort_inten);
    if isempty(maxInten)
        maxInten = 0;
    end
    
    tmp = sort_inten < 0.05 * maxInten; % Find results where intensity is less than 0.05 of the maximum abundance and discard them
    % sort_inten(tmp) = []; 
    sort_rts(tmp) = [];
    sort_ratioMatrix(tmp,:) = [];

    % Check if there are enough rows left
    % Replaces obj.hasMinRows(sort_ratioMatrix, minMSMSnum)
    if size(sort_ratioMatrix, 1) < minMSMSnum
        % If the ratio matrix has less than min rows, skip this group
        is_valid = false;
        sort_ratioMatrix = [];
    end
    % TODO: The filtering is not only here, but also the XIC integration area is determined later.
    %  Then, within the integration area, it is checked how many fall within this area.
    %  If there are not enough, they will be removed.
end
