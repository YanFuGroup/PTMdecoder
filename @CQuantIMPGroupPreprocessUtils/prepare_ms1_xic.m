function [rt_sorted, ratio_sorted, xic_rt, xic_intensity_smoothed, xic_intensity_raw, is_valid] = ...
    prepare_ms1_xic(cMs12DatasetIO, raw_name, rt_raw, intensity_raw, ratio_raw, ...
        minMSMSnum, low_mz_bound, high_mz_bound, selected_charge)
% Prepare MS1 inputs and load smoothed XIC.
% input:
%   cMs12DatasetIO (object)
%       MS12 dataset IO object
%   raw_name (1 x 1 char/string)
%       the name of the raw (mgf) file
%   rt_raw (N x 1 double) minutes
%       retention time in current group
%   intensity_raw (N x 1 double) intensity
%       intensity in current group
%   ratio_raw (N x K double)
%       ratio matrix of quantification in current group; rows aligned to rt_raw
%   minMSMSnum (1 x 1 double/int)
%       minimum MSMS number threshold
%   low_mz_bound (1 x 1 double) m/z
%       low precursor m/z bound
%   high_mz_bound (1 x 1 double) m/z
%       high precursor m/z bound
%   selected_charge (1 x 1 double/int)
%       current precursor charge
% output:
%   rt_sorted (N x 1 double) minutes
%       sorted retention times
%   ratio_sorted (N x K double)
%       sorted ratio matrix
%   xic_rt (M x 1 double) minutes
%       retention time grid
%   xic_intensity_smoothed (M x 1 double) intensity
%       smoothed intensity of XIC
%   xic_intensity_raw (M x 1 double) intensity
%       raw intensity of XIC
%   is_valid (1 x 1 logical)
%       is the input valid after preprocessing

[rt_sorted, ratio_sorted, is_valid] = ...
    CChromatogramUtils.preprocess_ms1_inputs(...
        rt_raw, intensity_raw, ratio_raw, minMSMSnum);

if ~is_valid
    xic_rt = [];
    xic_intensity_smoothed = [];
    xic_intensity_raw = [];
    return;
end

[xic_rt, xic_intensity_smoothed, xic_intensity_raw] = ...
    CChromatogramUtils.get_smoothed_xic(...
        cMs12DatasetIO, raw_name, low_mz_bound, high_mz_bound, selected_charge);
end
