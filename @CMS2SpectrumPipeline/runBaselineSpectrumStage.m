function [bSuccess,cstrIMP,abundance,ionTypePosCharge,ionIntens,frageff, ...
    is_X_not_full_column_rank,solver_diag,noise_model_fit_inputs,stability_cache] = ...
    runBaselineSpectrumStage(obj, peptideCtx, spectrumCtx)
% runBaselineSpectrumStage - Execute stage-1 baseline processing for one spectrum
% Input:
%   obj (CMS2SpectrumPipeline)
%       Processor instance carrying static solver/config parameters.
%   peptideCtx (1 x 1 struct)
%       Peptide-side context for one processing call.
%       Required fields:
%           pepSeq (1 x N char/string)
%               Peptide sequence.
%           isProtN (1 x 1 logical)
%               Whether peptide is protein N-terminal.
%           isProtC (1 x 1 logical)
%               Whether peptide is protein C-terminal.
%   spectrumCtx (1 x 1 struct)
%       Spectrum-side context for one processing call.
%       Required fields:
%           datasetName (1 x N char/string)
%               Dataset(raw file) name of the spectrum.
%           specName (1 x N char/string)
%               Spectrum title/index name.
%           expPeaks (P x 2 double)
%               Experimental peaks in [mz, intensity] form.
%           iCharge (1 x 1 double)
%               Precursor charge state.
%           precursorMZ (1 x 1 double)
%               Precursor m/z.
% Output:
%   bSuccess (1 x 1 logical)
%       True if this spectrum is processed successfully; false otherwise.
%   cstrIMP (Q x 1 cell)
%       Candidate peptidoform strings.
%   abundance (Q x 1 double)
%       Relative abundance for each peptidoform in cstrIMP.
%   ionTypePosCharge (R x 3 double)
%       Ion descriptors [ionType, ionPos, ionCharge].
%   ionIntens (R x 1 double)
%       Ion intensity corresponding to ionTypePosCharge.
%   frageff (R x 1 double)
%       Fragment efficiency corresponding to ionTypePosCharge.
%   is_X_not_full_column_rank (1 x 1 logical)
%       Whether matrix X is rank-deficient during abundance solve.
%   solver_diag (1 x 1 struct)
%       Stage-1 solver diagnostics placeholder. Fields:
%           is_X_not_full_column_rank (logical)
%               Whether matrix X is rank-deficient during abundance solve.
%           jaccard_stability (double)
%               Jaccard stability of the solution (NaN if not computed).
%           vif_all_imp_max (double)
%               Max VIF over all IMP columns in the solve design matrix.
%           vif_reported_imp_max (double)
%               Max VIF over filter-reported IMP columns.
%           vif_reported_each (double array)
%               Per-IMP VIF on filter-reported IMP columns.
%           support_frequency (double array)
%               Support frequency of each peptidoform (empty if not computed).
%           abundance_mad (double array)
%               Abundance MAD of each peptidoform (empty if not computed).
%           reported_imp_indices (double array)
%               Indices of peptidoforms reported as present (abundance >= tau).
%           num_successful_resamples (double)
%               Number of successful resamples if stability estimation is performed.
%   noise_model_fit_inputs (1 x 1 struct)
%       Dataset-level Noise Fitting Inputs.
%   stability_cache (1 x 1 struct)
%       Stage-1 spectrum cache placeholder for stage-2 stability estimation.

ionTypePosCharge = [];
ionIntens = [];
frageff = [];
is_X_not_full_column_rank = false;
solver_diag = struct( ...
    'is_X_not_full_column_rank', false, ...
    'jaccard_stability', NaN, ...
    'vif_all_imp_max', NaN, ...
    'vif_reported_imp_max', NaN, ...
    'vif_reported_each', [], ...
    'support_frequency', [], ...
    'abundance_mad', [], ...
    'reported_imp_indices', [], ...
    'num_successful_resamples', 0);
noise_model_fit_inputs = struct( ...
    'filteredOutExpPeakCount', 0, ...
    'filteredOutExpPeakSqSum', 0, ...
    'matchedExpPeaks', zeros(0, 3), ...
    'fittedMatchedPeakIntensities', zeros(0, 1));
stability_cache = struct( ...
    'vNonRedunTheoryIonMz', [], ...
    'matchedExpPeaks', zeros(0, 3), ...
    'massArrangement', [], ...
    'abundance', [], ...
    'fittedMatchedPeakIntensities', zeros(0, 1), ...
    'cstrIMP', {{}}, ...
    'solver_diag', solver_diag);

pepSeq = peptideCtx.pepSeq;
isProtN = peptideCtx.isProtN;
isProtC = peptideCtx.isProtC;
specName = spectrumCtx.specName;
expPeaks = spectrumCtx.expPeaks;
iCharge = spectrumCtx.iCharge;
dPrecursorMass = (spectrumCtx.precursorMZ - CConstant.pmass) * iCharge;
    
% Record dataset-level background summary from raw experimental peaks in
%   normalized-intensity space: low-intensity peaks below alpha are treated as
%   background candidates for sigma_base fitting.
if ~isempty(expPeaks)
    rawIntensity = expPeaks(:,2);
    rawIntensityNorm = rawIntensity / (max(rawIntensity) + eps);
    filteredOutMask = rawIntensityNorm < obj.m_alpha;
    filteredOutNormIntensity = rawIntensityNorm(filteredOutMask);
    noise_model_fit_inputs.filteredOutExpPeakCount = numel(filteredOutNormIntensity);
    noise_model_fit_inputs.filteredOutExpPeakSqSum = sum(filteredOutNormIntensity .^ 2);
end

% Build calculation context for mass-side modules
mass_ctx = struct( ...
    'm_pepSeq', pepSeq, ...
    'm_isProtN', isProtN, ...
    'm_isProtC', isProtC, ...
    'm_fixedModNameMass', {obj.m_fixedModNameMass}, ...
    'm_variableModNameMass', {obj.m_variableModNameMass}, ...
    'm_ms1_tolerance', obj.m_ms1_tolerance, ...
    'm_enzyme', obj.m_enzyme, ...
    'm_max_mod_per_peptide', obj.m_max_mod_per_peptide, ...
    'm_dPrecursorMass', dPrecursorMass, ...
    'm_strSpecName', specName, ...
    'm_iCharge', iCharge, ...
    'm_ionTypes', obj.m_ionTypes);

% Enumerate feasible peptidoforms and corresponding non-redundant ions
fixedPosMod = CMS2MassCalculator.getFixedPosMod(mass_ctx);
[bSuccess,inxSites,massArrangement] = CMS2MassCalculator.getMassArrangement(mass_ctx,fixedPosMod);
if ~bSuccess
    cstrIMP = [];
    abundance = [];
    return;
end

vNonRedunTheoryIonMz = CMS2MassCalculator.getNonRedunIons(mass_ctx,inxSites,massArrangement,fixedPosMod);

% Match and preprocess peaks on non-redundant ion space
matchedExpPeaks = CMS2PeakMatcher.match(expPeaks,vNonRedunTheoryIonMz,obj.m_ms2_tolerance);
% Apply alpha threshold to matched peaks for downstream processing
matchedExpPeaks = CMS2PeakMatcher.processPeaks(matchedExpPeaks,obj.m_alpha);

if isempty(matchedExpPeaks)
    if size(massArrangement,1) ~= 1
        cstrIMP = [];
        abundance = [];
        bSuccess = false;
        CLogger.debug(['[CMS2SpectrumPipeline:runBaselineSpectrumStage] ', ...
            'There is no non-redundant peak for discriminating peptidoforms for %s in %s.'], ...
            pepSeq, specName);
        return;
    else
        abundance = 1;
        cstrIMP = CMS2ResultIO.formatImpStrings(massArrangement,fixedPosMod,obj.m_variableModNameMass,inxSites,pepSeq);
        solver_diag.jaccard_stability = 1;
        solver_diag.vif_all_imp_max = 1;
        solver_diag.vif_reported_imp_max = 1;
        solver_diag.vif_reported_each = 1;
        solver_diag.reported_imp_indices = 1;
        solver_diag.support_frequency = 1;
        solver_diag.abundance_mad = 0;
        solver_diag.num_successful_resamples = 0;
        stability_cache.vNonRedunTheoryIonMz = vNonRedunTheoryIonMz;
        stability_cache.massArrangement = massArrangement;
        stability_cache.abundance = abundance;
        stability_cache.cstrIMP = cstrIMP;
        stability_cache.solver_diag = solver_diag;
        return;
    end
end

% Post-match peptidoform filtering (remove empty / strict-subset candidates)
[massArrangement, vNonRedunTheoryIonMz] = CMS2PeakMatcher.postMatchFilterPeptidoforms( ...
    matchedExpPeaks, massArrangement, vNonRedunTheoryIonMz);
if size(massArrangement,1) == 1
    abundance = 1;
    fittedMatchedPeakIntensities = CMS2QuantSolver.computeFittedMatchedPeakIntensities( ...
        vNonRedunTheoryIonMz, matchedExpPeaks, abundance);
    cstrIMP = CMS2ResultIO.formatImpStrings(massArrangement,fixedPosMod,obj.m_variableModNameMass,inxSites,pepSeq);
    bSuccess = true;
    solver_diag.jaccard_stability = 1;
    solver_diag.vif_all_imp_max = 1;
    solver_diag.vif_reported_imp_max = 1;
    solver_diag.vif_reported_each = 1;
    solver_diag.reported_imp_indices = 1;
    solver_diag.support_frequency = 1;
    solver_diag.abundance_mad = 0;
    solver_diag.num_successful_resamples = 0;
    noise_model_fit_inputs.matchedExpPeaks = matchedExpPeaks;
    noise_model_fit_inputs.fittedMatchedPeakIntensities = fittedMatchedPeakIntensities;
    stability_cache.vNonRedunTheoryIonMz = vNonRedunTheoryIonMz;
    stability_cache.matchedExpPeaks = matchedExpPeaks;
    stability_cache.massArrangement = massArrangement;
    stability_cache.abundance = abundance;
    stability_cache.fittedMatchedPeakIntensities = fittedMatchedPeakIntensities;
    stability_cache.cstrIMP = cstrIMP;
    stability_cache.solver_diag = solver_diag;
    return;
end

% Solve abundance according to configured model/method
solver_cfg = struct( ...
    'model', obj.m_model, ...
    'method', obj.m_method, ...
    'lambda', obj.m_lambda, ...
    'case_penalty_intens', obj.m_case_penalty_intens, ...
    'grid_penalty_intens', obj.m_grid_penalty_intens, ...
    'case_OLS_intens_weight', obj.m_case_OLS_intens_weight);

[abundance, frageff, ionTypePosCharge, ionIntens, is_X_not_full_column_rank, X] = ...
    CMS2QuantSolver.solve(vNonRedunTheoryIonMz, matchedExpPeaks, massArrangement, solver_cfg);

fittedMatchedPeakIntensities = CMS2QuantSolver.computeFittedMatchedPeakIntensities( ...
    vNonRedunTheoryIonMz, matchedExpPeaks, abundance, frageff, ionTypePosCharge);

% Final thresholding + normalization
reported_imp_mask = CMS2QuantSolver.getReportedImpMask(abundance, obj.m_resFilterThres);
abundance(~reported_imp_mask) = 0;
abundance = abundance / (sum(abundance) + eps);

cstrIMP = CMS2ResultIO.formatImpStrings(massArrangement,fixedPosMod,obj.m_variableModNameMass,inxSites,pepSeq);

solver_diag.is_X_not_full_column_rank = is_X_not_full_column_rank;
[solver_diag.vif_all_imp_max, solver_diag.vif_reported_imp_max, solver_diag.vif_reported_each] = CMS2QuantSolver.computeImpVifMetrics( ...
    X, size(massArrangement, 1), reported_imp_mask);
solver_diag.reported_imp_indices = find(reported_imp_mask);

noise_model_fit_inputs.matchedExpPeaks = matchedExpPeaks;
noise_model_fit_inputs.fittedMatchedPeakIntensities = fittedMatchedPeakIntensities;

stability_cache.vNonRedunTheoryIonMz = vNonRedunTheoryIonMz;
stability_cache.matchedExpPeaks = matchedExpPeaks;
stability_cache.massArrangement = massArrangement;
stability_cache.abundance = abundance;
stability_cache.fittedMatchedPeakIntensities = fittedMatchedPeakIntensities;
stability_cache.cstrIMP = cstrIMP;
stability_cache.solver_diag = solver_diag;

bSuccess = true;
end
