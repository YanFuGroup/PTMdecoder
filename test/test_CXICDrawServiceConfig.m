function tests = test_CXICDrawServiceConfig
% Validate CXICDrawServiceConfig stability filter parsing.
% Inputs:
%    (none)
% Outputs:
%    tests (matlab.unittest.Test)
%        Function-based test suite.

tests = functiontests(localfunctions);
end


function testDefaultStabilityFilterExists(testCase)
% Verify draw config exposes default MSMS stability filter.
% Inputs:
%    testCase (matlab.unittest.TestCase)
%        Test case context.
% Outputs:
%    (none)

param_map = makeBaseParamMap();
cfg = CXICDrawServiceConfig.fromParamMap(param_map);

testCase.verifyTrue(isfield(cfg, 'msms_stability_filter'));
testCase.verifyFalse(cfg.msms_stability_filter.enabled);
testCase.verifyTrue(cfg.msms_stability_filter.nan_as_fail);
testCase.verifyTrue(isfield(cfg, 'xic_layout'));
testCase.verifyEmpty(cfg.xic_layout);
end


function testOptionalXicLayoutPreserved(testCase)
% Verify draw config preserves an optional XIC layout override.
% Inputs:
%    testCase (matlab.unittest.TestCase)
%        Test case context.
% Outputs:
%    (none)

param_map = makeBaseParamMap();
xic_layout = struct();
xic_layout.figure_width_px = 1333;
xic_layout.figure_height_px = 667;
xic_layout.axes_width_fraction = 0.75;
param_map(CPTMdecoderWorkflowParamKeys.PARAM_XIC_LAYOUT) = xic_layout;

cfg = CXICDrawServiceConfig.fromParamMap(param_map);

testCase.verifyEqual(cfg.xic_layout.figure_width_px, 1333);
testCase.verifyEqual(cfg.xic_layout.figure_height_px, 667);
testCase.verifyEqual(cfg.xic_layout.axes_width_fraction, 0.75);
end


function param_map = makeBaseParamMap()
% Build a minimal valid parameter map for draw-XIC config.
% Inputs:
%    (none)
% Outputs:
%    param_map (containers.Map)
%        Valid key-value map for CXICDrawServiceConfig.

param_map = containers.Map('KeyType', 'char', 'ValueType', 'any');
param_map(CPTMdecoderWorkflowParamKeys.PARAM_MOD_FILE_PATH) = 'modify.ini';
param_map(CPTMdecoderWorkflowParamKeys.PARAM_FIXED_MOD) = '';
param_map(CPTMdecoderWorkflowParamKeys.PARAM_VARIABLE_MOD) = 'Acetyl[K]';
param_map(CPTMdecoderWorkflowParamKeys.PARAM_SPEC_DIR_PATH) = 'data/demo/spectra';
param_map(CPTMdecoderWorkflowParamKeys.PARAM_MS1_TOLERANCE_VALUE) = '10';
param_map(CPTMdecoderWorkflowParamKeys.PARAM_MS1_TOLERANCE_TYPE) = 'ppm';
param_map(CPTMdecoderWorkflowParamKeys.PARAM_ALPHA) = '0.01';
param_map(CPTMdecoderWorkflowParamKeys.PARAM_RESULT_FILTER_THRESHOLD) = '0.1';
param_map(CPTMdecoderWorkflowParamKeys.PARAM_OUTPUT_DIR_PATH) = 'output/demo';
param_map(CPTMdecoderWorkflowParamKeys.PARAM_CHECKED_PEPTIDES_RES_PATH) = 'output/demo/report_peptide_all_checked.txt';
param_map(CPTMdecoderWorkflowParamKeys.PARAM_MSMS_RES_PATH) = 'output/demo/report_msms.txt';
param_map(CPTMdecoderWorkflowParamKeys.PARAM_MIN_MSMS_NUM) = 1;
end
