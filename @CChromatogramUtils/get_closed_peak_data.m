function [rec_rt, rec_inten] = get_closed_peak_data(rt_grid, intensity_full, idx_start, idx_end)
% get_closed_peak_data
% Extracts XIC data with 1-point padding and zero-filled boundaries.
%
% Input:
%   rt_grid (N x 1 double) minutes
%       Full retention time grid
%   intensity_full (N x 1 double) intensity
%       Full intensity vector (column vector)
%   idx_start (1 x 1 double/int)
%       Original start index
%   idx_end (1 x 1 double/int)
%       Original end index
%
% Output:
%   rec_rt (M x 1 double) minutes
%       Reconstructed RT vector (padded)
%   rec_inten (M x 1 double) intensity
%       Reconstructed Intensity vector (padded with 0 at ends)

    % 1. Expand RT range if possible
    pad_start = max(1, idx_start - 1);
    pad_end   = min(length(rt_grid), idx_end + 1);
    
    rec_rt = rt_grid(pad_start : pad_end);
    
    % 2. Get Raw Intensity Segment
    % Get the actual intensities for the padded range first.
    rec_inten = intensity_full(pad_start : pad_end);
    
    % 3. Apply Zero-Padding (Clamping)
    % Only force zero if we actually expanded the range.
    % If we are at the boundary (e.g., idx_start=1), keep the real value.
    
    if pad_start < idx_start
        rec_inten(1) = 0; 
    end
    
    if pad_end > idx_end
        rec_inten(end) = 0;
    end
end
