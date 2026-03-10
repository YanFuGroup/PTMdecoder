function tests = test_CMS2MassCalculator
% TEST_CMS2MASSCALCULATOR Boundary tests for CMS2MassCalculator
% Input:
%   (none)
% Output:
%   tests (matlab.unittest.Test)

tests = functiontests(localfunctions);
end

function testTrypsinCTermBoundaryFilter(testCase)
% TESTTRYPSINCTERMBOUNDARYFILTER Validate trypsin C-term filtering removes invalid C-term masses
% Input:
%   testCase (matlab.unittest.TestCase)
% Output:
%   (none)

ctx = baseCtx('AK');
ctx.m_enzyme = struct('name', 'trypsin', 'limits', [1, 0]);
ctx.m_variableModNameMass = {
    'VarK2', 'K', 2.0
};
fixedPosMod = {};
neutralMass = CMS2MassCalculator.getNeutralPeptideTheoryMass(ctx, fixedPosMod);
ctx.m_dPrecursorMass = neutralMass + 2.0;
ctx.m_strSpecName = 'spec_trypsin_c_term_filter';

[bSuccess, ~, massArrangement] = CMS2MassCalculator.getMassArrangement(ctx, fixedPosMod);

testCase.verifyFalse(bSuccess);
testCase.verifyEmpty(massArrangement);
end

function testPpmVsDaToleranceBoundary(testCase)
% TESTPPMVSDATOLERANCEBOUNDARY Validate ppm/Da tolerance critical behavior on same mass mismatch
% Input:
%   testCase (matlab.unittest.TestCase)
% Output:
%   (none)

ctx = baseCtx('AA');
ctx.m_enzyme = struct('name', 'none', 'limits', []);
ctx.m_variableModNameMass = {
    'VarA1', 'A', 1.0
};
fixedPosMod = {};
neutralMass = CMS2MassCalculator.getNeutralPeptideTheoryMass(ctx, fixedPosMod);

% Keep a mismatch of 0.008 Da to probe tolerance boundary
ctx.m_dPrecursorMass = neutralMass + 1.008;
ctx.m_strSpecName = 'spec_tol_boundary';

ctx.m_ms1_tolerance = struct('value', 0.005, 'isppm', false);
[bDa, ~, massArrangementDa] = CMS2MassCalculator.getMassArrangement(ctx, fixedPosMod);

% Huge tolerance in ppm at this mass scale
ctx.m_ms1_tolerance = struct('value', 100, 'isppm', true);
[bPpm, ~, massArrangementPpm] = CMS2MassCalculator.getMassArrangement(ctx, fixedPosMod);

testCase.verifyFalse(bDa, 'Da tolerance should reject the 0.008 Da mismatch when tolerance=0.005 Da.');
testCase.verifyEmpty(massArrangementDa);

testCase.verifyTrue(bPpm, 'ppm tolerance should accept the mismatch under 10 ppm at precursor mass scale.');
testCase.verifyNotEmpty(massArrangementPpm);
end

function testFixedAndVariableCombination(testCase)
% TESTFIXEDANDVARIABLECOMBINATION Validate fixed/variable modification combination handling
% Input:
%   testCase (matlab.unittest.TestCase)
% Output:
%   (none)

ctx = baseCtx('AC');
ctx.m_isProtN = false;
ctx.m_isProtC = true;
ctx.m_enzyme = struct('name', 'none', 'limits', []);
ctx.m_fixedModNameMass = {
    'ProteinNMod', 'ProteinN-term', 10.0; ... % should be skipped because m_isProtN = false
    'CarbC', 'C', 57.02146; ...               % should be applied at C position
    'ProteinCMod', 'ProteinC-term', 3.0       % should be applied because m_isProtC = true
};
ctx.m_variableModNameMass = {
    'VarA1', 'A', 1.0
};

fixedPosMod = CMS2MassCalculator.getFixedPosMod(ctx);

testCase.verifyEqual(size(fixedPosMod, 1), 2, 'Expected two applicable fixed mods (C and ProteinC-term).');
testCase.verifyFalse(any(strcmp(fixedPosMod(:,2), 'ProteinNMod')));

theoryMassWithFixed = CMS2MassCalculator.getNeutralPeptideTheoryMass(ctx, fixedPosMod);
ctx.m_dPrecursorMass = theoryMassWithFixed + 1.0;
ctx.m_strSpecName = 'spec_fixed_var_comb';
ctx.m_ms1_tolerance = struct('value', 0.01, 'isppm', false);

[bSuccess, inxSites, massArrangement] = CMS2MassCalculator.getMassArrangement(ctx, fixedPosMod);

testCase.verifyTrue(bSuccess);
testCase.verifyNotEmpty(inxSites);
testCase.verifyNotEmpty(massArrangement);
testCase.verifyTrue(any(abs(massArrangement(:) - 1.0) < 1e-10), ...
    'Expected variable modification mass 1.0 to appear in arrangement.');
end

function ctx = baseCtx(pepSeq)
% BASECTX Build minimal context for CMS2MassCalculator tests

ctx = struct();
ctx.m_pepSeq = pepSeq;
ctx.m_isProtN = true;
ctx.m_isProtC = true;
ctx.m_fixedModNameMass = {};
ctx.m_variableModNameMass = {};
ctx.m_ms1_tolerance = struct('value', 0.01, 'isppm', false);
ctx.m_dPrecursorMass = 0;
ctx.m_enzyme = struct('name', 'none', 'limits', []);
ctx.m_strSpecName = 'spec';
ctx.m_iCharge = 2;
ctx.m_ionTypes = [1,2];
end
