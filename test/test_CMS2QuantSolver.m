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
