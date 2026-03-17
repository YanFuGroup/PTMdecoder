function perturbedMatchedExpPeaks = perturbMatchedPeaks(matchedExpPeaks, fittedMatchedPeakIntensities, noise_model, seed)
% perturbMatchedPeaks - Apply heteroscedastic perturbation to matched peaks
% Inputs:
%   matchedExpPeaks (K x 3 double)
%       Matched peaks with fixed column semantics:
%       col 1: ion index
%       col 2: normalized intensity (primary for stability perturbation)
%       col 3: compatibility raw-scale intensity.
%              The output col 3 is reconstructed from col 2 using
%              input-scale ratio max(col3)/(max(col2)+eps).
%   fittedMatchedPeakIntensities (K x 1 double)
%       Baseline fitted matched-peak intensities (yhat).
%   noise_model (struct)
%       Required fields:
%       - sigma_base (1 x 1 double): additive noise scale.
%       - gamma (1 x 1 double): multiplicative-relative noise scale.
%   seed (1 x 1 double, optional)
%       Random seed for deterministic perturbation.
% Outputs:
%   perturbedMatchedExpPeaks (K x 3 double)
%       Perturbed matched peaks preserving row count and ion index,
%       with col 3 restored to the input raw-intensity scale.

if nargin < 4
    seed = [];
end

if ~isnumeric(matchedExpPeaks) || size(matchedExpPeaks, 2) ~= 3
    CLogger.error(['[CMS2QuantSolver:InvalidMatchedExpPeaks] ', ...
        'matchedExpPeaks must be a numeric matrix with 3 columns.']);
end
if ~isnumeric(fittedMatchedPeakIntensities)
    CLogger.error(['[CMS2QuantSolver:InvalidFittedIntensities] ', ...
        'fittedMatchedPeakIntensities must be numeric.']);
end
if ~isstruct(noise_model) || ~isfield(noise_model, 'sigma_base') || ~isfield(noise_model, 'gamma')
    CLogger.error(['[CMS2QuantSolver:InvalidNoiseModel] ', ...
        'noise_model must contain sigma_base and gamma fields.']);
end
if ~isscalar(noise_model.sigma_base) || ~isnumeric(noise_model.sigma_base) || noise_model.sigma_base < 0
    CLogger.error(['[CMS2QuantSolver:InvalidSigmaBase] ', ...
        'noise_model.sigma_base must be a non-negative numeric scalar.']);
end
if ~isscalar(noise_model.gamma) || ~isnumeric(noise_model.gamma) || noise_model.gamma < 0
    CLogger.error(['[CMS2QuantSolver:InvalidGamma] ', ...
        'noise_model.gamma must be a non-negative numeric scalar.']);
end

fittedMatchedPeakIntensities = fittedMatchedPeakIntensities(:);
num_rows = size(matchedExpPeaks, 1);
if numel(fittedMatchedPeakIntensities) ~= num_rows
    CLogger.error(['[CMS2QuantSolver:SizeMismatch] ', ...
        'fittedMatchedPeakIntensities must have one value per matched peak row.']);
end

if ~isempty(seed)
    if ~isscalar(seed) || ~isnumeric(seed)
        CLogger.error('[CMS2QuantSolver:InvalidSeed] seed must be a numeric scalar.');
    end
    prev_rng = rng;
    cleanup = onCleanup(@() rng(prev_rng));
    rng(double(seed), 'twister');
end

perturbedMatchedExpPeaks = matchedExpPeaks;
if num_rows == 0
    return;
end

y = matchedExpPeaks(:, 2);
yhat = fittedMatchedPeakIntensities;
input_norm_max = max(y);
input_raw_max = max(matchedExpPeaks(:, 3));
raw_scale = input_raw_max / (input_norm_max + eps);
if ~isfinite(raw_scale) || raw_scale < 0
    raw_scale = 0;
end

% Compute heteroscedastic noise components:
% 1. Additive noise (eta): Baseline noise, scaled by sigma_base.
% 2. Multiplicative noise (epsilon): Relative to signal strength, scaled by gamma * yhat.
eta = noise_model.sigma_base * randn(num_rows, 1);
epsilon = noise_model.gamma * randn(num_rows, 1);
perturbed_y = max(y + eta + yhat .* epsilon, 0);

% Re-normalize intensities relative to the new maximum peak
scale = max(perturbed_y);
perturbed_norm = perturbed_y / (scale + eps);

perturbedMatchedExpPeaks(:, 2) = perturbed_norm;
perturbedMatchedExpPeaks(:, 3) = perturbed_norm * raw_scale;
end
