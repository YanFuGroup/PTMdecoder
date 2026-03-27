classdef CMS2QuantSolver
    % Static solver wrappers for MS2 quantification models

    methods (Static)
        noise_model = estimateDatasetNoiseModel(noise_model_fit_inputs)

        stability_diag = estimateStability(vNonRedunTheoryIonMz, matchedExpPeaks, massArrangement, solver_cfg, base_abundance, fittedMatchedPeakIntensities, noise_model, stability_options)

        [abundance, frageff, ionTypePosCharge, ionIntens, X] = solve(vNonRedunTheoryIonMz, matchedExpPeaks, massArrangement, solver_cfg)

        [vif_all_imp_max, vif_reported_imp_max, vif_reported_imp_each] = computeImpVifMetrics(X, num_imp, reported_imp_mask)

        reported_imp_mask = getReportedImpMask(abundance, relative_threshold)

        jaccard = computeJaccardIndex(maskA, maskB)

        perturbedMatchedExpPeaks = perturbMatchedPeaks(matchedExpPeaks, fittedMatchedPeakIntensities, noise_model, seed)

        fitted = computeFittedMatchedPeakIntensities(vNonRedunTheoryIonMz, matchedExpPeaks, abundance, frageff, ionTypePosCharge)

        penalty_factor = calculatePenaltyFactor(vNonRedunTheoryIonMz, matchedExpPeaks, lambda, case_penalty_intens, grid_penalty_intens)

        [X,ionTypePosCharge,ionIntens] = calculateX_FEV(vNonRedunTheoryIonMz, matchedExpPeaks, case_OLS_intens_weight)

        X = calculateX_Guan(vNonRedunTheoryIonMz)

        Y = calculateY_Guan(vNonRedunTheoryIonMz, matchedExpPeaks)

        Y = calculateY_1(vNonRedunTheoryIonMz, matchedExpPeaks)

        abundance = coreFEV_OLS(X, massArrangement)

        [abundance, frageff] = coreFEV_lasso(X, massArrangement, lambda)

        [abundance, frageff] = coreFEV_penalty(X, massArrangement, penalty_factor)

        abundance = coreGuan_OLS(X, Y)

        abundance = coreGuan_lasso(X, Y, lambda)

        abundance = coreGuan_penalty(X, Y, penalty_factor)
    end
end
