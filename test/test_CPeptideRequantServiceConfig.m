function tests = test_CPeptideRequantServiceConfig
% Validate CPeptideRequantServiceConfig stability filter parsing.
% Inputs:
%    (none)
% Outputs:
%    tests (matlab.unittest.Test)
%        Function-based test suite.

tests = functiontests(localfunctions);
end


function testDefaultStabilityFilterDisabled(testCase)
% Verify default MSMS stability filter config when params are omitted.
% Inputs:
%    testCase (matlab.unittest.TestCase)
%        Test case context.
% Outputs:
%    (none)

param_map = makeBaseParamMap();
cfg = CPeptideRequantServiceConfig.fromParamMap(param_map);

testCase.verifyTrue(isfield(cfg, 'msms_stability_filter'));
testCase.verifyFalse(cfg.msms_stability_filter.enabled);
testCase.verifyEmpty(cfg.msms_stability_filter.min_jaccard);
testCase.verifyEmpty(cfg.msms_stability_filter.min_support_frequency);
testCase.verifyEmpty(cfg.msms_stability_filter.max_abundance_mad);
testCase.verifyTrue(cfg.msms_stability_filter.nan_as_fail);
end


function testStabilityFilterParsesThresholds(testCase)
% Verify explicit stability thresholds are parsed correctly.
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

cfg = CPeptideRequantServiceConfig.fromParamMap(param_map);

testCase.verifyTrue(cfg.msms_stability_filter.enabled);
testCase.verifyEqual(cfg.msms_stability_filter.min_jaccard, 0.70, 'AbsTol', 1e-12);
testCase.verifyEqual(cfg.msms_stability_filter.min_support_frequency, 0.60, 'AbsTol', 1e-12);
testCase.verifyEqual(cfg.msms_stability_filter.max_abundance_mad, 0.15, 'AbsTol', 1e-12);
testCase.verifyFalse(cfg.msms_stability_filter.nan_as_fail);
end


function testInvalidMinJaccardThrowsError(testCase)
% Verify invalid jaccard threshold is rejected.
% Inputs:
%    testCase (matlab.unittest.TestCase)
%        Test case context.
% Outputs:
%    (none)

param_map = makeBaseParamMap();
param_map(CPTMdecoderWorkflowParamKeys.PARAM_PEPTIDE_QUANT_MIN_JACCARD_STABILITY) = '1.10';

f = @() CPeptideRequantServiceConfig.fromParamMap(param_map);
verifyLoggedErrorContains(testCase, f, 'CPeptideRequantServiceConfig:InvalidMinJaccardStability');
end


function testInvalidMaxMadThrowsError(testCase)
% Verify invalid abundance MAD threshold is rejected.
% Inputs:
%    testCase (matlab.unittest.TestCase)
%        Test case context.
% Outputs:
%    (none)

param_map = makeBaseParamMap();
param_map(CPTMdecoderWorkflowParamKeys.PARAM_PEPTIDE_QUANT_MAX_ABUNDANCE_MAD) = '-0.1';

f = @() CPeptideRequantServiceConfig.fromParamMap(param_map);
verifyLoggedErrorContains(testCase, f, 'CPeptideRequantServiceConfig:InvalidMaxAbundanceMad');
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
% Build a minimal valid parameter map for peptide requant config.
% Inputs:
%    (none)
% Outputs:
%    param_map (containers.Map)
%        Valid key-value map for CPeptideRequantServiceConfig.

param_map = containers.Map('KeyType', 'char', 'ValueType', 'any');
param_map(CPTMdecoderWorkflowParamKeys.PARAM_MOD_FILE_PATH) = 'modify.ini';
param_map(CPTMdecoderWorkflowParamKeys.PARAM_FIXED_MOD) = '';
param_map(CPTMdecoderWorkflowParamKeys.PARAM_VARIABLE_MOD) = 'Acetyl[K]';
param_map(CPTMdecoderWorkflowParamKeys.PARAM_SPEC_DIR_PATH) = 'data/demo/spectra';
param_map(CPTMdecoderWorkflowParamKeys.PARAM_MS1_TOLERANCE_VALUE) = '10';
param_map(CPTMdecoderWorkflowParamKeys.PARAM_MS1_TOLERANCE_TYPE) = 'ppm';
param_map(CPTMdecoderWorkflowParamKeys.PARAM_ALPHA) = '0.01';
param_map(CPTMdecoderWorkflowParamKeys.PARAM_FASTA_FILE_PATH) = 'data/demo/database/demo.fasta';
param_map(CPTMdecoderWorkflowParamKeys.PARAM_REGULAR_EXPRESS) = '>([^ ,]*)';
param_map(CPTMdecoderWorkflowParamKeys.PARAM_FILTERED_RES_FILE_PATH) = 'data/demo/filtered_result.txt';
param_map(CPTMdecoderWorkflowParamKeys.PARAM_RESULT_FILTER_THRESHOLD) = '0.1';
param_map(CPTMdecoderWorkflowParamKeys.PARAM_OUTPUT_DIR_PATH) = 'output/demo';
param_map(CPTMdecoderWorkflowParamKeys.PARAM_CHECKED_PEPTIDES_RES_PATH) = 'output/demo/report_peptide_all_checked.txt';
param_map(CPTMdecoderWorkflowParamKeys.PARAM_MSMS_RES_PATH) = 'output/demo/report_msms.txt';
end
