function [vif_all_imp_max, vif_reported_imp_max, vif_reported_imp_each] = computeImpVifMetrics(X, num_imp, reported_imp_mask)
% Compute max VIF on all IMP columns and reported IMP columns
% Inputs:
%   X (N x P double)
%       Solver design matrix.
%   num_imp (1 x 1 double)
%       Number of IMP columns at the beginning of X.
%   reported_imp_mask (num_imp x 1 logical)
%       Reported IMP mask after abundance filtering.
% Outputs:
%   vif_all_imp_max (1 x 1 double)
%       Max VIF over all IMP columns.
%   vif_reported_imp_max (1 x 1 double)
%       Max VIF over reported IMP columns.
%   vif_reported_imp_each (R x 1 double)
%       Per-IMP VIF vector on reported IMP columns.

if nargin < 3
    CLogger.error(['[CMS2QuantSolver:MissingVifInputs] ', ...
        'X, num_imp, and reported_imp_mask must be provided.']);
end

if isempty(X) || ~isnumeric(X) || ~ismatrix(X)
    CLogger.error('[CMS2QuantSolver:InvalidVifDesignMatrix] X must be a non-empty numeric 2-D matrix.');
end
if size(X, 2) == 0
    CLogger.error('[CMS2QuantSolver:InvalidVifDesignMatrix] X must contain at least one column.');
end

if ~isscalar(num_imp) || ~isnumeric(num_imp) || num_imp <= 0 || mod(num_imp, 1) ~= 0
    CLogger.error('[CMS2QuantSolver:InvalidNumImp] num_imp must be a positive integer scalar.');
end
if num_imp > size(X, 2)
    CLogger.error(['[CMS2QuantSolver:InvalidNumImp] ', ...
        'num_imp (%d) cannot exceed size(X,2) (%d).'], num_imp, size(X, 2));
end

if isempty(reported_imp_mask)
    CLogger.error('[CMS2QuantSolver:MissingReportedImpMask] reported_imp_mask must be provided.');
end
reported_imp_mask = logical(reported_imp_mask(:));
if numel(reported_imp_mask) ~= num_imp
    CLogger.error('[CMS2QuantSolver:InvalidReportedImpMaskLength] reported_imp_mask length must equal num_imp.');
end

X_imp = X(:, 1:num_imp);
vif_all_imp_each = computeVifVector(X_imp);
vif_all_imp_max = max(vif_all_imp_each);

if nnz(reported_imp_mask) == 0
    vif_reported_imp_max = NaN;
    vif_reported_imp_each = NaN(0, 1);
else
    X_reported = X_imp(:, reported_imp_mask);
    vif_reported_imp_each = computeVifVector(X_reported);
    vif_reported_imp_max = max(vif_reported_imp_each);
end
end


function vif_vec = computeVifVector(X_sub)
% Compute uncentered VIF vector for a feature matrix (no-intercept model).
% X_sub should contain only the free variables (non-zero features) from the QP solution.

num_cols = size(X_sub, 2);

% 1. Handle the case where no features are selected (all penalized to 0)
if num_cols == 0
    vif_vec = NaN(0, 1);
    return;
end

% 2. Handle the case with only 1 free variable (no collinearity is possible)
if num_cols == 1
    vif_vec = 1;
    return;
end

% 3. Fast computation of uncentered VIF
Q = X_sub' * X_sub;

% Prevent execution crash from exact multicollinearity (singular matrix).
% rcond() estimates the reciprocal condition number; near 0 indicates rank deficiency.
if rcond(Q) < 1e-12 
    vif_vec = Inf(num_cols, 1);
else
    % Core mathematical equivalence: uncentered VIF
    % Derived from block matrix inversion (Schur complement)
    vif_vec = diag(Q) .* diag(inv(Q));
end

end