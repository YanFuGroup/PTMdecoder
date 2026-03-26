function [abundance, frageff, ionTypePosCharge, ionIntens, is_X_not_full_column_rank, X] = solve(vNonRedunTheoryIonMz, matchedExpPeaks, massArrangement, solver_cfg)
% solve - Model/method dispatch for MS2 quantification
% Inputs:
%   vNonRedunTheoryIonMz (L x T double)
%       Non-redundant ion table.
%   matchedExpPeaks (K x 3 double)
%       Matched/processed peaks [ion_index, normalized_intensity, raw_intensity].
%   massArrangement (M x S double)
%       Candidate peptidoform mass arrangements.
%   solver_cfg (struct)
%       Required fields: model, method, lambda, case_penalty_intens,
%       grid_penalty_intens, case_OLS_intens_weight.
% Outputs:
%   abundance (M x 1 double)
%       Relative abundance for each peptidoform.
%   frageff (Q x 1 double)
%       Fragmentation efficiency coefficients (FEV only; empty for FEC/FEE).
%   ionTypePosCharge (U x 3 double)
%       Ion tuple [type, position, charge] participating in X matrix.
%   ionIntens (U x 1 double)
%       Aggregated ion intensities corresponding to ionTypePosCharge.
%   is_X_not_full_column_rank (1 x 1 logical)
%       True if X is rank-deficient.
%   X (N x P double)
%       Design matrix used by the abundance solver.

ionTypePosCharge = [];
ionIntens = [];
frageff = [];
is_X_not_full_column_rank = false;
X = [];

if solver_cfg.method == 3
    % Penalty: compute peptidoform-dependent penalty factors.
    penalty_factor = CMS2QuantSolver.calculatePenaltyFactor(vNonRedunTheoryIonMz, matchedExpPeaks, solver_cfg.lambda, ...
        solver_cfg.case_penalty_intens, solver_cfg.grid_penalty_intens);
end

switch solver_cfg.model
    case 1
        % FEV model: variable fragmentation efficiency.
        [X,ionTypePosCharge,ionIntens] = CMS2QuantSolver.calculateX_FEV(vNonRedunTheoryIonMz, matchedExpPeaks, ...
            solver_cfg.case_OLS_intens_weight);
        switch solver_cfg.method
            case 1
                abundance = CMS2QuantSolver.coreFEV_OLS(X, massArrangement);
            case 2
                [abundance,frageff] = CMS2QuantSolver.coreFEV_lasso(X, massArrangement, solver_cfg.lambda);
            case 3
                [abundance,frageff] = CMS2QuantSolver.coreFEV_penalty(X, massArrangement, penalty_factor);
        end
    case 2
        % FEC model: constant-but-unequal fragmentation efficiency.
        X = CMS2QuantSolver.calculateX_Guan(vNonRedunTheoryIonMz);
        Y = CMS2QuantSolver.calculateY_Guan(vNonRedunTheoryIonMz, matchedExpPeaks);
        switch solver_cfg.method
            case 1
                abundance = CMS2QuantSolver.coreGuan_OLS(X, Y);
            case 2
                abundance = CMS2QuantSolver.coreGuan_lasso(X, Y, solver_cfg.lambda);
            case 3
                abundance = CMS2QuantSolver.coreGuan_penalty(X, Y, penalty_factor);
        end
    case 3
        % FEE model: equal fragmentation efficiency.
        X = CMS2QuantSolver.calculateX_Guan(vNonRedunTheoryIonMz);
        Y = CMS2QuantSolver.calculateY_1(vNonRedunTheoryIonMz, matchedExpPeaks);
        switch solver_cfg.method
            case 1
                abundance = CMS2QuantSolver.coreGuan_OLS(X, Y);
            case 2
                abundance = CMS2QuantSolver.coreGuan_lasso(X, Y, solver_cfg.lambda);
        end
end

if rank(X)~=size(X,2)
    is_X_not_full_column_rank=true;
end
end
