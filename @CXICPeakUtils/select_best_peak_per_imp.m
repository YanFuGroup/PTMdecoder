function idx_selected = select_best_peak_per_imp(imp_max_props, area_imp_by_peak)
% select_best_peak_per_imp
% Select the best XIC peak for each IMP based on a weighted score.
%
% Inputs:
%   imp_max_props (K x P double)
%       Max ratio contribution per IMP per peak
%   area_imp_by_peak (K x P double) area
%       Area contribution per IMP per peak
%
% Output:
%   idx_selected (K x 1 double)
%       Best peak index per IMP

num_imp = size(imp_max_props, 1);
idx_selected = zeros(num_imp, 1);

for idx_imp = 1:num_imp
    % Get features for this IMP across all peaks
    % Ratio (Purity)
    props = imp_max_props(idx_imp, :);
    % Area Contribution
    areas = area_imp_by_peak(idx_imp, :);

    % Calculate Weighted Score = Area * Ratio
    % This handles two cases:
    % 1. Co-eluting isomers (70:30): Small ratio difference, huge area difference -> Select Main Peak.
    % 2. Distinct elution (1% in Big Peak vs 95% in Small Peak):
    %    Big Peak Score: High Area * Very Low Ratio -> Low Score
    %    Small Peak Score: Low Area * High Ratio -> High Score (due to ratio squared effect effectively) -> Select Small Peak.
    scores = areas .* props;
    [~, best_idx] = max(scores);
    idx_selected(idx_imp) = best_idx;
end
end
