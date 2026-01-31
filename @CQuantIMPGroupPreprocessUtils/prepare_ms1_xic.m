function [sort_rts, sort_ratioMatrix, rt_grid, smoothed_intensity, intensity, is_valid] = ...
    prepare_ms1_xic(cMs12DatasetIO, raw_name, current_rts, current_inten, current_ratioMatrix, ...
        minMSMSnum, low_mz_bound, high_mz_bound, selected_charge)
% Prepare MS1 inputs and load smoothed XIC.
% input:
%   cMs12DatasetIO (object)
%       MS12 dataset IO object
%   raw_name (1 x 1 char/string)
%       the name of the raw (mgf) file
%   current_rts (N x 1 double) minutes
%       retention time in current group
%   current_inten (N x 1 double) intensity
%       intensity in current group
%   current_ratioMatrix (N x K double)
%       ratio matrix of quantification in current group; rows aligned to current_rts
%   minMSMSnum (1 x 1 double/int)
%       minimum MSMS number threshold
%   low_mz_bound (1 x 1 double) m/z
%       low precursor m/z bound
%   high_mz_bound (1 x 1 double) m/z
%       high precursor m/z bound
%   selected_charge (1 x 1 double/int)
%       current precursor charge
% output:
%   sort_rts (N x 1 double) minutes
%       sorted retention times
%   sort_ratioMatrix (N x K double)
%       sorted ratio matrix
%   rt_grid (M x 1 double) minutes
%       retention time grid
%   smoothed_intensity (M x 1 double) intensity
%       smoothed intensity of XIC
%   intensity (M x 1 double) intensity
%       raw intensity of XIC
%   is_valid (1 x 1 logical)
%       is the input valid after preprocessing

[sort_rts, sort_ratioMatrix, is_valid] = ...
    CChromatogramUtils.preprocess_ms1_inputs(...
        current_rts, current_inten, current_ratioMatrix, minMSMSnum);

if ~is_valid
    rt_grid = [];
    smoothed_intensity = [];
    intensity = [];
    return;
end

[rt_grid, smoothed_intensity, intensity] = ...
    CChromatogramUtils.get_smoothed_xic(...
        cMs12DatasetIO, raw_name, low_mz_bound, high_mz_bound, selected_charge);
end
