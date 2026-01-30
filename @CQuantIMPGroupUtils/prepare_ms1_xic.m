function [sort_rts, sort_ratioMatrix, rt_grid, smoothed_intensity, intensity, is_valid] = ...
    prepare_ms1_xic(cMs12DatasetIO, raw_name, current_rts, current_inten, current_ratioMatrix, ...
        minMSMSnum, low_mz_bound, high_mz_bound, selected_charge)
% Prepare MS1 inputs and load smoothed XIC.
% input:
%   cMs12DatasetIO
%       MS12 dataset IO object
%   raw_name
%       the name of the raw (mgf) file
%   current_rts
%       retention time in current group
%   current_inten
%       intensity in current group
%   current_ratioMatrix
%       ratio matrix of quantification in current group
%   minMSMSnum
%       minimum MSMS number threshold
%   low_mz_bound
%       low precursor m/z bound
%   high_mz_bound
%       high precursor m/z bound
%   selected_charge
%       current precursor charge
% output:
%   sort_rts
%       sorted retention times
%   sort_ratioMatrix
%       sorted ratio matrix
%   rt_grid
%       retention time grid
%   smoothed_intensity
%       smoothed intensity of XIC
%   intensity
%       raw intensity of XIC
%   is_valid
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
