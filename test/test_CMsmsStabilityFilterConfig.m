function tests = test_CMsmsStabilityFilterConfig
% Validate CMsmsStabilityFilterConfig parsing and validation.
% Inputs:
%    (none)
% Outputs:
%    tests (matlab.unittest.Test)
%        Function-based test suite.

tests = functiontests(localfunctions);
end


function testDefaultFilterConfig(testCase)
% Verify default filter config values when params are omitted.
% Inputs:
%    testCase (matlab.unittest.TestCase)
%        Test case context.
% Outputs:
%    (none)

param_map = makeBaseParamMap();
filter_cfg = CMsmsStabilityFilterConfig.fromParamMap(param_map, 'TestScope');

testCase.verifyFalse(filter_cfg.enabled);
testCase.verifyEmpty(filter_cfg.min_jaccard);
testCase.verifyEmpty(filter_cfg.min_support_frequency);
testCase.verifyEmpty(filter_cfg.max_abundance_mad);
testCase.verifyTrue(filter_cfg.nan_as_fail);
end


function testExplicitFilterConfig(testCase)
% Verify explicit filter params are parsed correctly.
% Inputs:
%    testCase (matlab.unittest.TestCase)
%        Test case context.
% Outputs:
%    (none)

param_map = makeBaseParamMap();
param_map(CPTMdecoderWorkflowParamKeys.PARAM_PEPTIDE_QUANT_STABILITY_FILTER_ON) = '1';
param_map(CPTMdecoderWorkflowParamKeys.PARAM_PEPTIDE_QUANT_MIN_JACCARD_STABILITY) = '0.70';
param_map(CPTMdecoderWorkflowParamKeys.PARAM_PEPTIDE_QUANT_MIN_SUPPORT_FREQUENCY) = '0.60';
param_map(CPTMdecoderWorkflowParamKeys.PARAM_PEPTIDE_QUANT_MAX_ABUNDANCE_MAD) = '0.15';
param_map(CPTMdecoderWorkflowParamKeys.PARAM_PEPTIDE_QUANT_STABILITY_NAN_AS_FAIL) = '0';

filter_cfg = CMsmsStabilityFilterConfig.fromParamMap(param_map, 'TestScope');

testCase.verifyTrue(filter_cfg.enabled);
testCase.verifyEqual(filter_cfg.min_jaccard, 0.70, 'AbsTol', 1e-12);
testCase.verifyEqual(filter_cfg.min_support_frequency, 0.60, 'AbsTol', 1e-12);
testCase.verifyEqual(filter_cfg.max_abundance_mad, 0.15, 'AbsTol', 1e-12);
testCase.verifyFalse(filter_cfg.nan_as_fail);
end


function testInvalidMinJaccardThrowsScopedError(testCase)
% Verify min jaccard range validation uses caller-provided error prefix.
% Inputs:
%    testCase (matlab.unittest.TestCase)
%        Test case context.
% Outputs:
%    (none)

param_map = makeBaseParamMap();
param_map(CPTMdecoderWorkflowParamKeys.PARAM_PEPTIDE_QUANT_MIN_JACCARD_STABILITY) = '1.10';

f = @() CMsmsStabilityFilterConfig.fromParamMap(param_map, 'MyCaller');
verifyLoggedErrorContains(testCase, f, 'MyCaller:InvalidMinJaccardStability');
end


function testInvalidMaxMadThrowsScopedError(testCase)
% Verify max abundance MAD validation uses caller-provided error prefix.
% Inputs:
%    testCase (matlab.unittest.TestCase)
%        Test case context.
% Outputs:
%    (none)

param_map = makeBaseParamMap();
param_map(CPTMdecoderWorkflowParamKeys.PARAM_PEPTIDE_QUANT_MAX_ABUNDANCE_MAD) = '-0.1';

f = @() CMsmsStabilityFilterConfig.fromParamMap(param_map, 'MyCaller');
verifyLoggedErrorContains(testCase, f, 'MyCaller:InvalidMaxAbundanceMad');
end


function verifyLoggedErrorContains(testCase, funcHandle, expectedTag)
% Verify CLogger error id and expected business tag.
% Inputs:
%    testCase (matlab.unittest.TestCase)
%        Test case context.
%    funcHandle (1 x 1 function_handle)
%        Function expected to fail through CLogger.error.
%    expectedTag (1 x N char/string)
%        Business error tag in message.
% Outputs:
%    (none)

did_throw = false;
caught_error = [];
try
    funcHandle();
catch me
    did_throw = true;
    caught_error = me;
end
testCase.assertTrue(did_throw, 'Expected function to raise an error.');
testCase.assertEqual(caught_error.identifier, 'CLogger:LoggedError');
testCase.verifyTrue(contains(caught_error.message, expectedTag), ...
    ['Expected error message to contain tag: ', expectedTag]);
end


function param_map = makeBaseParamMap()
% Build base param map for stability filter parser tests.
% Inputs:
%    (none)
% Outputs:
%    param_map (containers.Map)
%        Parameter map used by parser tests.

param_map = containers.Map('KeyType', 'char', 'ValueType', 'any');
end
