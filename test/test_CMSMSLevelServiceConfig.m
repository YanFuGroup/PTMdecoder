function tests = test_CMSMSLevelServiceConfig
% TEST_CMSMSLEVELSERVICECONFIG Unit tests for CMSMSLevelServiceConfig.fromParamMap.
% Input:
%   (none)
% Output:
%   tests (matlab.unittest.Test)

tests = functiontests(localfunctions);
end


function testMissingNResamplesSkipsStabilityOptions(testCase)
% Verify missing n_resamples keeps Stage-2 disabled at config layer.
% Inputs:
%    testCase (matlab.unittest.TestCase)
%        Test case context.
% Outputs:
%    (none)

param_map = makeBaseParamMap();
cfg = CMSMSLevelServiceConfig.fromParamMap(param_map);

testCase.verifyFalse(isfield(cfg, 'stability_options'));
end


function testOnlyNResamplesUsesDefaultRandomSeed(testCase)
% Verify n_resamples alone enables Stage-2 with default random_seed=42.
% Inputs:
%    testCase (matlab.unittest.TestCase)
%        Test case context.
% Outputs:
%    (none)

param_map = makeBaseParamMap();
param_map(CPTMdecoderWorkflowParamKeys.PARAM_N_RESAMPLES) = '30';

cfg = CMSMSLevelServiceConfig.fromParamMap(param_map);

testCase.verifyTrue(isfield(cfg, 'stability_options'));
testCase.verifyEqual(cfg.stability_options.n_resamples, 30);
testCase.verifyEqual(cfg.stability_options.random_seed, 42);
testCase.verifyEqual(cfg.stability_options.relative_threshold, cfg.result_filter_threshold);
end


function testInvalidNResamplesThrowsError(testCase)
% Verify non-positive or non-integer n_resamples is rejected.
% Inputs:
%    testCase (matlab.unittest.TestCase)
%        Test case context.
% Outputs:
%    (none)

param_map = makeBaseParamMap();
param_map(CPTMdecoderWorkflowParamKeys.PARAM_N_RESAMPLES) = '2.5';

f = @() CMSMSLevelServiceConfig.fromParamMap(param_map);
testCase.verifyError(f, 'CMSMSLevelServiceConfig:InvalidNResamples');
end


function testRandomSeedWithoutNResamplesDoesNotEnableStage2(testCase)
% Verify random_seed alone does not enable Stage-2 config.
% Inputs:
%    testCase (matlab.unittest.TestCase)
%        Test case context.
% Outputs:
%    (none)

param_map = makeBaseParamMap();
param_map(CPTMdecoderWorkflowParamKeys.PARAM_RANDOM_SEED) = '7';

cfg = CMSMSLevelServiceConfig.fromParamMap(param_map);

testCase.verifyFalse(isfield(cfg, 'stability_options'));
end


function testExplicitRandomSeedIsRespected(testCase)
% Verify explicit random_seed is preserved when Stage-2 is enabled.
% Inputs:
%    testCase (matlab.unittest.TestCase)
%        Test case context.
% Outputs:
%    (none)

param_map = makeBaseParamMap();
param_map(CPTMdecoderWorkflowParamKeys.PARAM_N_RESAMPLES) = '12';
param_map(CPTMdecoderWorkflowParamKeys.PARAM_RANDOM_SEED) = '7';

cfg = CMSMSLevelServiceConfig.fromParamMap(param_map);

testCase.verifyTrue(isfield(cfg, 'stability_options'));
testCase.verifyEqual(cfg.stability_options.n_resamples, 12);
testCase.verifyEqual(cfg.stability_options.random_seed, 7);
testCase.verifyEqual(cfg.stability_options.relative_threshold, cfg.result_filter_threshold);
end


function param_map = makeBaseParamMap()
% Build a minimal valid parameter map for CMSMSLevelServiceConfig.
% Inputs:
%    (none)
% Outputs:
%    param_map (containers.Map)
%        Valid key-value map for fromParamMap.

param_map = containers.Map('KeyType', 'char', 'ValueType', 'any');
param_map(CPTMdecoderWorkflowParamKeys.PARAM_MOD_FILE_PATH) = 'modify.ini';
param_map(CPTMdecoderWorkflowParamKeys.PARAM_FIXED_MOD) = '';
param_map(CPTMdecoderWorkflowParamKeys.PARAM_VARIABLE_MOD) = 'Acetyl[K]';
param_map(CPTMdecoderWorkflowParamKeys.PARAM_SPEC_DIR_PATH) = 'data/demo/spectra';
param_map(CPTMdecoderWorkflowParamKeys.PARAM_MS1_TOLERANCE_VALUE) = '10';
param_map(CPTMdecoderWorkflowParamKeys.PARAM_MS1_TOLERANCE_TYPE) = 'ppm';
param_map(CPTMdecoderWorkflowParamKeys.PARAM_MS2_TOLERANCE) = '0.02';
param_map(CPTMdecoderWorkflowParamKeys.PARAM_ALPHA) = '0.01';
param_map(CPTMdecoderWorkflowParamKeys.PARAM_FASTA_FILE_PATH) = 'data/demo/database/demo.fasta';
param_map(CPTMdecoderWorkflowParamKeys.PARAM_REGULAR_EXPRESS) = '>([^ ,]*)';
param_map(CPTMdecoderWorkflowParamKeys.PARAM_FILTERED_RES_FILE_PATH) = 'data/demo/filtered_result.txt';
param_map(CPTMdecoderWorkflowParamKeys.PARAM_MODEL) = '1';
param_map(CPTMdecoderWorkflowParamKeys.PARAM_METHOD) = '1';
param_map(CPTMdecoderWorkflowParamKeys.PARAM_RESULT_FILTER_THRESHOLD) = '0.1';
param_map(CPTMdecoderWorkflowParamKeys.PARAM_ENZYME_NAME) = 'trypsin';
param_map(CPTMdecoderWorkflowParamKeys.PARAM_ENZYME_LIMIT_C_TERM_POSSIBLE_MOD) = '14.015650';
param_map(CPTMdecoderWorkflowParamKeys.PARAM_OUTPUT_DIR_PATH) = 'output/demo';
param_map(CPTMdecoderWorkflowParamKeys.PARAM_PEP_SPEC_FILE_PATH) = 'data/demo/pepSpecFile.txt';
end
