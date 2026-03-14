classdef CMS2QuantSolver
    % Static solver wrappers for MS2 quantification models

    methods (Static)
        [abundance, frageff, ionTypePosCharge, ionIntens, is_X_not_full_column_rank] = solve(vNonRedunTheoryIonMz, matchedExpPeaks, massArrangement, solver_cfg)

        [reported_imp_mask, tau] = getReportedImpMask(abundance, relative_threshold)

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
