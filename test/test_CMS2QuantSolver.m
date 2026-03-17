function tests = test_CMS2QuantSolver
% TEST_CMS2QUANTSOLVER Test script for CMS2QuantSolver solve dispatch behavior
% Input:
%   (none)
% Output:
%   tests (matlab.unittest.Test)

tests = functiontests(localfunctions);
end


function testSolveModel2Method1Basic(testCase)
% TESTSOLVEMODEL2METHOD1BASIC Validate model=2, method=1 returns normalized abundance and rank flag false
% Input:
%   testCase (matlab.unittest.TestCase)
% Output:
%   (none)

testCase.assumeNotEmpty(which('quadprog'), 'Optimization Toolbox (quadprog) is required for this test.');

v = [
    100, 1, 1, 1, 0, 1, 1, 0;
    200, 2, 2, 1, 0, 2, 0, 1
];
matched = [
    1, 0.6, 60;
    2, 0.4, 40
];
massArrangement = [0; 1];
solver_cfg = struct('model',2,'method',1,'lambda',0.1, ...
    'case_penalty_intens','intens_sum','grid_penalty_intens','intens_sum','case_OLS_intens_weight','none');

[abundance, frageff, ionTypePosCharge, ionIntens, isRankDef] = CMS2QuantSolver.solve(v, matched, massArrangement, solver_cfg);

testCase.verifySize(abundance, [2,1]);
testCase.verifyEqual(sum(abundance), 1, 'AbsTol', 1e-6);
testCase.verifyEmpty(frageff);
testCase.verifyEmpty(ionTypePosCharge);
testCase.verifyEmpty(ionIntens);
testCase.verifyFalse(isRankDef);
end


function testSolveRankDeficientX(testCase)
% TESTSOLVERANKDEFICIENTX Validate rank-deficient X is reported
% Input:
%   testCase (matlab.unittest.TestCase)
% Output:
%   (none)

testCase.assumeNotEmpty(which('quadprog'), 'Optimization Toolbox (quadprog) is required for this test.');

v = [
    100, 1, 1, 1, 0, 1, 1, 1;
    200, 2, 2, 1, 0, 2, 1, 1
];
matched = [
    1, 0.5, 50;
    2, 0.5, 50
];
massArrangement = [0; 1];
solver_cfg = struct('model',2,'method',1,'lambda',0.1, ...
    'case_penalty_intens','intens_sum','grid_penalty_intens','intens_sum','case_OLS_intens_weight','none');

[abundance, ~, ~, ~, isRankDef] = CMS2QuantSolver.solve(v, matched, massArrangement, solver_cfg);

testCase.verifySize(abundance, [2,1]);
testCase.verifyTrue(isRankDef);
end


function testGetReportedImpMaskUsesRelativeThreshold(testCase)
% Validate reported IMP mask uses tau=0.01*max(abundance)
% Input:
%   testCase (matlab.unittest.TestCase)
% Output:
%   (none)

abundance = [0; 0.02; 1.0];
reportedMask = CMS2QuantSolver.getReportedImpMask(abundance, 0.01);

testCase.verifyEqual(reportedMask, logical([false; true; true]));
end


function testGetReportedImpMaskIncludesBoundaryAtTau(testCase)
% Validate abundance exactly equal to tau is included (>= tau)
% Input:
%   testCase (matlab.unittest.TestCase)
% Output:
%   (none)

abundance = [0.01; 0.5; 1.0];
reportedMask = CMS2QuantSolver.getReportedImpMask(abundance, 0.01);

testCase.verifyEqual(reportedMask, logical([true; true; true]));
end


function testGetReportedImpMaskExclude(testCase)
% Validate abundance less than tau is excluded
% Input:
%   testCase (matlab.unittest.TestCase)
% Output:
%   (none)

abundance = [0.01; 0.5; 1.0];
reportedMask = CMS2QuantSolver.getReportedImpMask(abundance, 0.1);

testCase.verifyEqual(reportedMask, logical([false; true; true]));
end


function testGetReportedImpMaskWithZeroThresholdKeepsNonNegative(testCase)
% Validate zero threshold keeps all entries with abundance >= 0
% Input:
%   testCase (matlab.unittest.TestCase)
% Output:
%   (none)

abundance = [0; 0.1; 0.3];
reportedMask = CMS2QuantSolver.getReportedImpMask(abundance, 0);

testCase.verifyEqual(reportedMask, logical([true; true; true]));
end


function testGetReportedImpMaskAllZeroReturnsAllFalse(testCase)
% Validate all-zero abundance returns all false without throwing
% Input:
%   testCase (matlab.unittest.TestCase)
% Output:
%   (none)

abundance = [0; 0; 0];

reportedMask = CMS2QuantSolver.getReportedImpMask(abundance, 0.01);

testCase.verifyEqual(reportedMask, false(3, 1));
end


function testGetReportedImpMaskAllNegativeReturnsAllFalse(testCase)
% Validate all-negative abundance returns all false without throwing
% Input:
%   testCase (matlab.unittest.TestCase)
% Output:
%   (none)

abundance = [-1; -0.1; -2];

reportedMask = CMS2QuantSolver.getReportedImpMask(abundance, 0.01);

testCase.verifyEqual(reportedMask, false(3, 1));
end


function testComputeJaccardIndexBasicCases(testCase)
% Validate identical/disjoint/empty-mask behaviors
% Input:
%   testCase (matlab.unittest.TestCase)
% Output:
%   (none)

maskA = logical([true; false; true; false]);
maskB = logical([true; false; true; false]);
maskC = logical([false; true; false; true]);
emptyA = false(5,1);
emptyB = false(5,1);

testCase.verifyEqual(CMS2QuantSolver.computeJaccardIndex(maskA, maskB), 1, 'AbsTol', 1e-12);
testCase.verifyEqual(CMS2QuantSolver.computeJaccardIndex(maskA, maskC), 0, 'AbsTol', 1e-12);
testCase.verifyEqual(CMS2QuantSolver.computeJaccardIndex(emptyA, emptyB), 0, 'AbsTol', 1e-12);
end


function testPerturbMatchedPeaksInvariants(testCase)
% Validate shape and value invariants after perturbation
% Input:
%   testCase (matlab.unittest.TestCase)
% Output:
%   (none)

matched = [
    3, 0.30, 30;
    5, 0.60, 60;
    9, 1.00, 100
];
yhat = [28; 58; 95];
noise_model = struct('sigma_base', 2.0, 'gamma', 0.1);

perturbed = CMS2QuantSolver.perturbMatchedPeaks(matched, yhat, noise_model, 1);

testCase.verifySize(perturbed, size(matched));
testCase.verifyEqual(perturbed(:,1), matched(:,1));
testCase.verifyGreaterThanOrEqual(perturbed(:,3), zeros(size(perturbed,1),1));
testCase.verifyGreaterThanOrEqual(perturbed(:,2), zeros(size(perturbed,1),1));
testCase.verifyLessThanOrEqual(perturbed(:,2), ones(size(perturbed,1),1));
expectedScale = max(matched(:,3)) / (max(matched(:,2)) + eps);
testCase.verifyEqual(perturbed(:,3), perturbed(:,2) * expectedScale, 'AbsTol', 1e-12);

if any(perturbed(:,3) > 0)
    testCase.verifyEqual(max(perturbed(:,2)), 1, 'AbsTol', 1e-12);
end
end


function testPerturbMatchedPeaksSeedReproducible(testCase)
% Validate fixed seed gives deterministic perturbation
% Input:
%   testCase (matlab.unittest.TestCase)
% Output:
%   (none)

matched = [
    1, 0.20, 20;
    2, 0.40, 40;
    3, 0.80, 80
];
yhat = [18; 39; 82];
noise_model = struct('sigma_base', 1.0, 'gamma', 0.05);

perturbed1 = CMS2QuantSolver.perturbMatchedPeaks(matched, yhat, noise_model, 1234);
perturbed2 = CMS2QuantSolver.perturbMatchedPeaks(matched, yhat, noise_model, 1234);

testCase.verifyEqual(perturbed1, perturbed2, 'AbsTol', 1e-12);
end


function testComputeFittedMatchedPeakIntensitiesPrimaryPath(testCase)
% Validate primary path uses site-specific fragmentation efficiency
% Input:
%   testCase (matlab.unittest.TestCase)
% Output:
%   (none)

v = [
    100, 1, 1, 1, 0, 0, 1, 0;
    200, 2, 2, 1, 0, 0, 0, 1;
    300, 1, 2, 1, 0, 0, 1, 1
];
matched = [
    1, 0.5, 50;
    2, 0.4, 40
];
abundance = [0.6; 0.4];
frageff = [0.8; 0.5];
ionTypePosCharge = [
    1, 1, 1;
    2, 2, 1
];

fitted = CMS2QuantSolver.computeFittedMatchedPeakIntensities(v, matched, abundance, frageff, ionTypePosCharge);

testCase.verifySize(fitted, [2, 1]);
testCase.verifyGreaterThanOrEqual(fitted, zeros(2, 1));
% fittedShape = [0.6*0.8; 0.4*0.5] = [0.48; 0.2]
% expIntensities = [0.5; 0.4] (from matched(:,2))
% c = (0.48*0.5 + 0.2*0.4) / (0.48^2 + 0.2^2) = 0.32 / 0.2704 = 1.18343195
testCase.verifyEqual(fitted(1), 0.48 * (0.32 / 0.2704), 'AbsTol', 1e-6);
testCase.verifyEqual(fitted(2), 0.20 * (0.32 / 0.2704), 'AbsTol', 1e-6);
end

function testComputeFittedMatchedPeakIntensitiesFallbackPath(testCase)
% Validate fallback path uses global scaling regression when frageff is missing
% Input:
%   testCase (matlab.unittest.TestCase)
% Output:
%   (none)

v = [
    100, 1, 1, 1, 0, 0, 1, 0;
    200, 2, 2, 1, 0, 0, 0, 1
];
matched = [
    1, 0.8, 80;
    2, 0.2, 20
];
abundance = [0.5; 0.5];

fitted = CMS2QuantSolver.computeFittedMatchedPeakIntensities(v, matched, abundance);

testCase.verifySize(fitted, [2, 1]);
testCase.verifyGreaterThanOrEqual(fitted, zeros(2, 1));
end

function testComputeFittedMatchedPeakIntensitiesEmpty(testCase)
% Validate empty input returns empty output safely
% Input:
%   testCase (matlab.unittest.TestCase)
% Output:
%   (none)

v = [100, 1, 1, 1, 0, 0, 1, 0];
matched = [];
abundance = [1.0; 0.0];

fitted = CMS2QuantSolver.computeFittedMatchedPeakIntensities(v, matched, abundance);

testCase.verifyEmpty(fitted);
end


function testEstimateStabilityDeterministicPositiveSupportCase(testCase)
% Validate deterministic zero-noise case where reported IMP support is strictly positive
% Input:
%   testCase (matlab.unittest.TestCase)
% Output:
%   (none)

testCase.assumeNotEmpty(which('quadprog'), 'Optimization Toolbox (quadprog) is required for this test.');

v = [
    100, 1, 1, 1, 0, 1, 1, 0, 0;
    200, 2, 2, 1, 0, 2, 0, 1, 0
];
matched = [
    1, 0.6, 60;
    2, 0.4, 40
];
massArrangement = [0; 1; 2];
solver_cfg = struct('model',2,'method',1,'lambda',0.1, ...
    'case_penalty_intens','intens_sum','grid_penalty_intens','intens_sum','case_OLS_intens_weight','none');

% Third IMP is structural zero in baseline and should never be reported.
base_abundance = [0.7; 0.3; 0];
fitted = [0.6; 0.4];
% Use zero perturbation noise to keep all successful resamples deterministic.
noise_model = struct('sigma_base', 0.0, 'gamma', 0.0, 'tau_floor', 0.1);
stability_options = struct('n_resamples', 5, 'random_seed', 1, 'relative_threshold', 0.01);

stability_diag = CMS2QuantSolver.estimateStability( ...
    v, matched, massArrangement, solver_cfg, base_abundance, fitted, noise_model, stability_options);

testCase.verifyGreaterThanOrEqual(stability_diag.jaccard_stability, 0);
testCase.verifyLessThanOrEqual(stability_diag.jaccard_stability, 1);
testCase.verifyEqual(stability_diag.num_successful_resamples, 5);
testCase.verifyEqual(stability_diag.reported_imp_indices, [1; 2]);
testCase.verifyFalse(any(stability_diag.reported_imp_indices == 3));
testCase.verifySize(stability_diag.support_frequency, [2, 1]);
testCase.verifyGreaterThan(stability_diag.support_frequency, zeros(2, 1));
testCase.verifyLessThanOrEqual(stability_diag.support_frequency, ones(2, 1));
testCase.verifySize(stability_diag.abundance_mad, [2, 1]);
testCase.verifyGreaterThanOrEqual(stability_diag.abundance_mad, zeros(2, 1));
end


function testEstimateStabilityNonZeroNoiseSupportCanBeZero(testCase)
% Validate that with non-zero noise, a baseline-reported IMP can have zero support
% Input:
%   testCase (matlab.unittest.TestCase)
% Output:
%   (none)

testCase.assumeNotEmpty(which('quadprog'), 'Optimization Toolbox (quadprog) is required for this test.');

% IMP #2 is structurally unsupported by X (all zeros in IMP-2 column),
% so it can be baseline-reported but still have zero support in resamples.
v = [
    100, 1, 1, 1, 0, 1, 1, 0;
    200, 2, 2, 1, 0, 2, 1, 0
];
matched = [
    1, 0.6, 60;
    2, 0.4, 40
];
massArrangement = [0; 1];
solver_cfg = struct('model',2,'method',1,'lambda',0.1, ...
    'case_penalty_intens','intens_sum','grid_penalty_intens','intens_sum','case_OLS_intens_weight','none');

base_abundance = [0.9; 0.1];
fitted = [0.6; 0.4];
noise_model = struct('sigma_base', 0.02, 'gamma', 0.05, 'tau_floor', 0.1);
stability_options = struct('n_resamples', 6, 'random_seed', 7, 'relative_threshold', 0.01);

stability_diag = CMS2QuantSolver.estimateStability( ...
    v, matched, massArrangement, solver_cfg, base_abundance, fitted, noise_model, stability_options);

testCase.verifyEqual(stability_diag.reported_imp_indices, [1; 2]);
testCase.verifySize(stability_diag.support_frequency, [2, 1]);
testCase.verifyGreaterThanOrEqual(stability_diag.support_frequency, zeros(2, 1));
testCase.verifyLessThanOrEqual(stability_diag.support_frequency, ones(2, 1));
testCase.verifyEqual(stability_diag.support_frequency(2), 0, 'AbsTol', 1e-12);
testCase.verifySize(stability_diag.abundance_mad, [2, 1]);
testCase.verifyGreaterThanOrEqual(stability_diag.abundance_mad, zeros(2, 1));
end


function testEstimateStabilitySeedReproducible(testCase)
% Validate estimateStability is reproducible with fixed random_seed
% Input:
%   testCase (matlab.unittest.TestCase)
% Output:
%   (none)

testCase.assumeNotEmpty(which('quadprog'), 'Optimization Toolbox (quadprog) is required for this test.');

v = [
    100, 1, 1, 1, 0, 1, 1, 0;
    200, 2, 2, 1, 0, 2, 0, 1
];
matched = [
    1, 0.6, 60;
    2, 0.4, 40
];
massArrangement = [0; 1];
solver_cfg = struct('model',2,'method',1,'lambda',0.1, ...
    'case_penalty_intens','intens_sum','grid_penalty_intens','intens_sum','case_OLS_intens_weight','none');

base_abundance = [0.7; 0.3];
fitted = [0.6; 0.4];
noise_model = struct('sigma_base', 0.02, 'gamma', 0.05, 'tau_floor', 0.1);
stability_options = struct('n_resamples', 6, 'random_seed', 11, 'relative_threshold', 0.01);

diag1 = CMS2QuantSolver.estimateStability( ...
    v, matched, massArrangement, solver_cfg, base_abundance, fitted, noise_model, stability_options);
diag2 = CMS2QuantSolver.estimateStability( ...
    v, matched, massArrangement, solver_cfg, base_abundance, fitted, noise_model, stability_options);

testCase.verifyEqual(diag1.jaccard_stability, diag2.jaccard_stability, 'AbsTol', 1e-12);
testCase.verifyEqual(diag1.num_successful_resamples, diag2.num_successful_resamples);
testCase.verifyEqual(diag1.reported_imp_indices, diag2.reported_imp_indices);
testCase.verifyEqual(diag1.support_frequency, diag2.support_frequency, 'AbsTol', 1e-12);
testCase.verifyEqual(diag1.abundance_mad, diag2.abundance_mad, 'AbsTol', 1e-12);
end


function testEstimateStabilityZeroNoiseHasPerfectStability(testCase)
% Validate zero-noise perturb-and-resolve yields perfect stability diagnostics
% Inputs:
%   testCase (matlab.unittest.TestCase)
%       MATLAB unit test handle.
% Outputs:
%   (none)

testCase.assumeNotEmpty(which('quadprog'), 'Optimization Toolbox (quadprog) is required for this test.');

v = [
    100, 1, 1, 1, 0, 1, 1, 0;
    200, 2, 2, 1, 0, 2, 0, 1
];
matched = [
    1, 0.6, 60;
    2, 0.4, 40
];
massArrangement = [0; 1];
solver_cfg = struct('model',2,'method',1,'lambda',0.1, ...
    'case_penalty_intens','intens_sum','grid_penalty_intens','intens_sum','case_OLS_intens_weight','none');

base_abundance = [0.7; 0.3];
fitted = [0.6; 0.4];
noise_model = struct('sigma_base', 0.0, 'gamma', 0.0, 'tau_floor', 0.1);
stability_options = struct('n_resamples', 7, 'random_seed', 13, 'relative_threshold', 0.01);

stability_diag = CMS2QuantSolver.estimateStability( ...
    v, matched, massArrangement, solver_cfg, base_abundance, fitted, noise_model, stability_options);

testCase.verifyEqual(stability_diag.num_successful_resamples, stability_options.n_resamples);
testCase.verifyEqual(stability_diag.reported_imp_indices, [1; 2]);
testCase.verifyEqual(stability_diag.jaccard_stability, 1, 'AbsTol', 1e-12);
testCase.verifyEqual(stability_diag.support_frequency, ones(2, 1), 'AbsTol', 1e-12);
testCase.verifyEqual(stability_diag.abundance_mad, zeros(2, 1), 'AbsTol', 1e-12);
end


function testEstimateStabilityTooManyFailuresThrows(testCase)
% Validate estimateStability raises when failed resamples exceed half
% Input:
%   testCase (matlab.unittest.TestCase)
% Output:
%   (none)

testCase.assumeNotEmpty(which('quadprog'), 'Optimization Toolbox (quadprog) is required for this test.');

v = [
    100, 1, 1, 1, 0, 1, 1, 0;
    200, 2, 2, 1, 0, 2, 0, 1
];
matched = [
    1, 0.6, 60;
    2, 0.4, 40
];
massArrangement = [0; 1];
solver_cfg_invalid = struct('model',99,'method',1,'lambda',0.1, ...
    'case_penalty_intens','intens_sum','grid_penalty_intens','intens_sum','case_OLS_intens_weight','none');

base_abundance = [0.7; 0.3];
fitted = [0.6; 0.4];
noise_model = struct('sigma_base', 0.02, 'gamma', 0.05, 'tau_floor', 0.1);
stability_options = struct('n_resamples', 5, 'random_seed', 3, 'relative_threshold', 0.01);

f = @() CMS2QuantSolver.estimateStability( ...
    v, matched, massArrangement, solver_cfg_invalid, base_abundance, fitted, noise_model, stability_options);

testCase.verifyError(f, 'CLogger:LoggedError');
try
    f();
catch ME
    testCase.verifyThat(string(ME.message), matlab.unittest.constraints.ContainsSubstring('CMS2QuantSolver:TooManyResampleFailures'));
end
end


function testEstimateStabilityEmptyBaselineReportedSetThrows(testCase)
% Validate baseline empty reported IMP set is a hard failure
% Input:
%   testCase (matlab.unittest.TestCase)
% Output:
%   (none)

testCase.assumeNotEmpty(which('quadprog'), 'Optimization Toolbox (quadprog) is required for this test.');

v = [
    100, 1, 1, 1, 0, 1, 1, 0;
    200, 2, 2, 1, 0, 2, 0, 1
];
matched = [
    1, 0.6, 60;
    2, 0.4, 40
];
massArrangement = [0; 1];
solver_cfg = struct('model',2,'method',1,'lambda',0.1, ...
    'case_penalty_intens','intens_sum','grid_penalty_intens','intens_sum','case_OLS_intens_weight','none');

base_abundance = [0; 0];
fitted = [0.6; 0.4];
noise_model = struct('sigma_base', 0.02, 'gamma', 0.05, 'tau_floor', 0.1);
stability_options = struct('n_resamples', 5, 'random_seed', 3, 'relative_threshold', 0.01);

f = @() CMS2QuantSolver.estimateStability( ...
    v, matched, massArrangement, solver_cfg, base_abundance, fitted, noise_model, stability_options);

testCase.verifyError(f, 'CLogger:LoggedError');
try
    f();
catch ME
    testCase.verifyThat(string(ME.message), matlab.unittest.constraints.ContainsSubstring('CMS2QuantSolver:InvalidBaselineReportedIMP'));
end
end

