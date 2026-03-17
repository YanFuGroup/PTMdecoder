function stability_diag = estimateStability(vNonRedunTheoryIonMz, matchedExpPeaks, massArrangement, solver_cfg, base_abundance, fittedMatchedPeakIntensities, noise_model, stability_options)
% estimateStability - Execute perturb-and-resolve stability resampling
% Inputs:
%   vNonRedunTheoryIonMz (L x T double)
%       Non-redundant ion table used by the abundance solver.
%   matchedExpPeaks (K x 3 double)
%       Matched peaks with fixed column semantics:
%       col 1: ion index
%       col 2: normalized intensity
%       col 3: compatibility raw-scale column restored from col 2
%              using the input max(col3)/max(col2) ratio.
%   massArrangement (M x S double)
%       Candidate peptidoform mass arrangements. Size M must match length of base_abundance.
%   solver_cfg (struct)
%       Solver configuration for the baseline solve path.
%       Required fields:
%       - model (1 x 1 double): solver model id.
%       - method (1 x 1 double): solver method id.
%       - lambda (1 x 1 double): regularization parameter.
%       - case_penalty_intens (char/string): penalty strategy.
%       - grid_penalty_intens (char/string): penalty grid strategy.
%       - case_OLS_intens_weight (char/string): OLS weighting strategy.
%   base_abundance (M x 1 double)
%       Baseline abundance vector after stage-1 processing.
%   fittedMatchedPeakIntensities (K x 1 double)
%       Baseline fitted matched-peak intensities (yhat).
%   noise_model (struct)
%       Dataset-level noise model.
%       Required fields:
%       - sigma_base (1 x 1 double): additive noise scale.
%       - gamma (1 x 1 double): multiplicative-relative noise scale.
%   stability_options (struct)
%       Stability resampling options (must be provided by caller).
%       Required fields:
%       - n_resamples (1 x 1 double): number of resamples (positive integer).
%       - random_seed (1 x 1 double): base random seed.
%       - relative_threshold (1 x 1 double): reported-IMP threshold factor.
% Outputs:
%   stability_diag (struct)
%       Stability diagnostics.
%       Fields:
%       - jaccard_stability (1 x 1 double): mean Jaccard over successful resamples.
%       - support_frequency (R x 1 double): support frequency for baseline reported IMPs.
%       - abundance_mad (R x 1 double): abundance MAD over successful resamples.
%       - reported_imp_indices (R x 1 double): baseline reported IMP indices.
%       - num_successful_resamples (1 x 1 double): successful resample count.

% Parse strict caller-provided options (no local fallback defaults).
[num_resamples, base_seed, relative_threshold] = parseStabilityOptions(stability_options);

% Validate core dimensional consistency before entering the resampling loop.
if size(massArrangement, 1) ~= numel(base_abundance)
    CLogger.error('[CMS2QuantSolver:BaseAbundanceSizeMismatch] base_abundance length must match size(massArrangement,1).');
end
if size(matchedExpPeaks, 1) ~= numel(fittedMatchedPeakIntensities)
    CLogger.error('[CMS2QuantSolver:FittedSizeMismatch] fittedMatchedPeakIntensities length must match matchedExpPeaks rows.');
end

% Build baseline reported IMP set and initialize accumulators.
[base_reported_mask, ~] = CMS2QuantSolver.getReportedImpMask(base_abundance, relative_threshold);
reported_imp_indices = find(base_reported_mask);
if isempty(reported_imp_indices)
    CLogger.error(['[CMS2QuantSolver:InvalidBaselineReportedIMP] ', ...
        'Baseline reported IMP set is empty; skip stability estimation for this spectrum.']);
end
support_counts = zeros(numel(reported_imp_indices), 1);
resampled_abundance_reported = zeros(numel(reported_imp_indices), num_resamples);

sum_jaccard = 0;
num_successful_resamples = 0;
num_failed_resamples = 0;
max_failed_resamples_threshold = floor(num_resamples / 2);

% Perturb-and-resolve loop with failure-tolerant accounting.
for idxResample = 1:num_resamples
    resample_failed = false;
    is_exception_failure = false;
    exception_identifier = '';
    exception_message = '';
    try
        perturbedMatchedExpPeaks = CMS2QuantSolver.perturbMatchedPeaks( ...
            matchedExpPeaks, fittedMatchedPeakIntensities, noise_model, base_seed + idxResample);
        abundance_resampled = CMS2QuantSolver.solve( ...
            vNonRedunTheoryIonMz, perturbedMatchedExpPeaks, massArrangement, solver_cfg);
        [resampled_reported_mask, ~] = CMS2QuantSolver.getReportedImpMask(abundance_resampled, relative_threshold);

        if ~any(resampled_reported_mask)
            resample_failed = true;
        else
            sum_jaccard = sum_jaccard + CMS2QuantSolver.computeJaccardIndex(base_reported_mask, resampled_reported_mask);
            if ~isempty(reported_imp_indices)
                support_counts = support_counts + double(resampled_reported_mask(reported_imp_indices));
                resampled_abundance_reported(:, num_successful_resamples + 1) = abundance_resampled(reported_imp_indices);
            end
            num_successful_resamples = num_successful_resamples + 1;
        end
    catch ME
        resample_failed = true;
        is_exception_failure = true;
        exception_identifier = ME.identifier;
        exception_message = ME.message;
    end

    if resample_failed
        num_failed_resamples = num_failed_resamples + 1;
        if is_exception_failure
            CLogger.debug(['[CMS2QuantSolver:ResampleFailed] ', ...
                'resample %d/%d failed: [%s] %s'], ...
                idxResample, num_resamples, exception_identifier, exception_message);
        end
        if num_failed_resamples > max_failed_resamples_threshold
            CLogger.error(['[CMS2QuantSolver:TooManyResampleFailures] ', ...
                'Resample failures exceeded half of total runs (%d/%d).'], ...
                num_failed_resamples, num_resamples);
        end
    end
end

% Finalize metrics from successful runs only.
if num_successful_resamples > 0
    jaccard_stability = sum_jaccard / num_successful_resamples;
    support_frequency = support_counts / num_successful_resamples;
    abundance_mad = computeAbundanceMad(resampled_abundance_reported(:, 1:num_successful_resamples));
else
    jaccard_stability = NaN;
    support_frequency = NaN(size(support_counts));
    abundance_mad = NaN(size(support_counts));
end

stability_diag = struct( ...
    'jaccard_stability', jaccard_stability, ...
    'support_frequency', support_frequency, ...
    'abundance_mad', abundance_mad, ...
    'reported_imp_indices', reported_imp_indices, ...
    'num_successful_resamples', num_successful_resamples);

CLogger.debug(['[CMS2QuantSolver:StabilitySummary] B=%d, success=%d, failed=%d, ', ...
    'jaccard=%.6f, reported_count=%d'], ...
    num_resamples, num_successful_resamples, num_failed_resamples, ...
    stability_diag.jaccard_stability, numel(reported_imp_indices));
end


function [num_resamples, base_seed, relative_threshold] = parseStabilityOptions(stability_options)
% parseStabilityOptions - Validate strict stability options from caller
% Inputs:
%   stability_options (struct)
%       Stability options provided by caller.
%       Required fields:
%       - n_resamples (1 x 1 double): number of resamples (positive integer).
%       - random_seed (1 x 1 double): base seed for deterministic resampling.
%       - relative_threshold (1 x 1 double): reported-IMP threshold factor.
% Outputs:
%   num_resamples (1 x 1 double)
%       Number of resamples.
%   base_seed (1 x 1 double)
%       Base random seed.
%   relative_threshold (1 x 1 double)
%       Relative threshold used by reported-IMP mask extraction.

if ~isstruct(stability_options)
    CLogger.error('[CMS2QuantSolver:InvalidStabilityOptions] stability_options must be a struct.');
end
if ~isfield(stability_options, 'n_resamples')
    CLogger.error('[CMS2QuantSolver:MissingNResamples] stability_options.n_resamples is required.');
end
num_resamples = stability_options.n_resamples;
if ~isscalar(num_resamples) || ~isnumeric(num_resamples) || num_resamples <= 0 || mod(num_resamples, 1) ~= 0
    CLogger.error('[CMS2QuantSolver:InvalidNumResamples] n_resamples must be a positive integer.');
end
if ~isfield(stability_options, 'relative_threshold')
    CLogger.error('[CMS2QuantSolver:MissingRelativeThreshold] stability_options.relative_threshold is required.');
end
relative_threshold = stability_options.relative_threshold;
if ~isscalar(relative_threshold) || ~isnumeric(relative_threshold) || relative_threshold < 0
    CLogger.error('[CMS2QuantSolver:InvalidRelativeThreshold] relative_threshold must be a non-negative scalar.');
end
if ~isfield(stability_options, 'random_seed')
    CLogger.error('[CMS2QuantSolver:MissingRandomSeed] stability_options.random_seed is required.');
end
base_seed = stability_options.random_seed;
if ~isscalar(base_seed) || ~isnumeric(base_seed)
    CLogger.error('[CMS2QuantSolver:InvalidRandomSeed] random_seed must be a numeric scalar.');
end
end


function abundance_mad = computeAbundanceMad(resampled_abundance_reported)
% computeAbundanceMad - Compute per-IMP MAD on successful resample abundances
% Inputs:
%   resampled_abundance_reported (R x B_success double)
%       Resampled abundance matrix on baseline reported IMP index space.
% Outputs:
%   abundance_mad (R x 1 double)
%       Median absolute deviation for each baseline reported IMP.

num_imp = size(resampled_abundance_reported, 1);
abundance_mad = zeros(num_imp, 1);
for idxImp = 1:num_imp
    abundance_vec = resampled_abundance_reported(idxImp, :);
    med_val = median(abundance_vec);
    abundance_mad(idxImp) = median(abs(abundance_vec - med_val));
end
end