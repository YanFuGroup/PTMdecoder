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
[reportedMask, tau] = CMS2QuantSolver.getReportedImpMask(abundance, 0.01);

testCase.verifyEqual(tau, 0.01, 'AbsTol', 1e-12);
testCase.verifyEqual(reportedMask, logical([false; true; true]));
end


function testGetReportedImpMaskRejectsAllZero(testCase)
% Validate all-zero abundance is rejected
% Input:
%   testCase (matlab.unittest.TestCase)
% Output:
%   (none)

abundance = [0; 0; 0];

f = @() CMS2QuantSolver.getReportedImpMask(abundance, 0.01);
testCase.verifyError(f, 'CLogger:LoggedError');
try
    f();
catch ME
    testCase.verifyThat(string(ME.message), matlab.unittest.constraints.ContainsSubstring('CMS2QuantSolver:InvalidMaxAbundance'));
end
end


function testGetReportedImpMaskRejectsAllNegative(testCase)
% Validate all-negative abundance is rejected
% Input:
%   testCase (matlab.unittest.TestCase)
% Output:
%   (none)

abundance = [-1; -0.1; -2];

f = @() CMS2QuantSolver.getReportedImpMask(abundance, 0.01);
testCase.verifyError(f, 'CLogger:LoggedError');
try
    f();
catch ME
    testCase.verifyThat(string(ME.message), matlab.unittest.constraints.ContainsSubstring('CMS2QuantSolver:InvalidMaxAbundance'));
end
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
