function auxic = get_final_area_from_peak_areas(area_each_XIC_peak, idx_selected)
% get_final_area_from_peak_areas
% Get final area per IMP from precomputed peak areas.
%
% Inputs:
%   area_each_XIC_peak (K x P double) area
%       Area contribution per IMP per peak
%   idx_selected (K x 1 double)
%       Selected peak index per IMP
%
% Output:
%   auxic (K x 1 double) area
%       Final area per IMP

num_imp = size(area_each_XIC_peak, 1);
auxic = zeros(num_imp, 1);

for idx_imp = 1:num_imp
    sel_peak_idx = idx_selected(idx_imp);
    auxic(idx_imp, 1) = area_each_XIC_peak(idx_imp, sel_peak_idx);
end
end
