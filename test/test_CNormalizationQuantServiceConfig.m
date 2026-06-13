function tests = test_CNormalizationQuantServiceConfig
% Validate normalization quant sparse-XIC threshold parsing.
tests = functiontests(localfunctions);
end

function testMinXicNonzeroPointsDefaultsAndParses(testCase)
[param_map, cleanup] = makeBaseParamMap(); %#ok<NASGU>
cfg = CNormalizationQuantServiceConfig.fromParamMap(param_map);
testCase.verifyEqual(cfg.msms_cfg.min_xic_nonzero_points, 5);

param_map(CPTMdecoderWorkflowParamKeys.PARAM_MIN_XIC_NONZERO_POINTS) = '3';
cfg = CNormalizationQuantServiceConfig.fromParamMap(param_map);
testCase.verifyEqual(cfg.msms_cfg.min_xic_nonzero_points, 3);
end

function testMinXicNonzeroPointsRejectsInvalidValues(testCase)
[param_map, cleanup] = makeBaseParamMap(); %#ok<NASGU>
invalid_values = {'0', '-1', '2.5'};
for idx = 1:numel(invalid_values)
    param_map(CPTMdecoderWorkflowParamKeys.PARAM_MIN_XIC_NONZERO_POINTS) = invalid_values{idx};
    verifyLoggedErrorContains(testCase, ...
        @() CNormalizationQuantServiceConfig.fromParamMap(param_map), ...
        'CNormalizationQuantServiceConfig:InvalidMinXicNonzeroPoints');
end
end

function testMinMsmsNumDefaultsAndParses(testCase)
[param_map, cleanup] = makeBaseParamMap(); %#ok<NASGU>
cfg = CNormalizationQuantServiceConfig.fromParamMap(param_map);
testCase.verifyEqual(cfg.msms_cfg.min_MSMS_num, 1);

param_map(CPTMdecoderWorkflowParamKeys.PARAM_MIN_MSMS_NUM) = '';
cfg = CNormalizationQuantServiceConfig.fromParamMap(param_map);
testCase.verifyEqual(cfg.msms_cfg.min_MSMS_num, 1);

param_map(CPTMdecoderWorkflowParamKeys.PARAM_MIN_MSMS_NUM) = '3';
cfg = CNormalizationQuantServiceConfig.fromParamMap(param_map);
testCase.verifyEqual(cfg.msms_cfg.min_MSMS_num, 3);
end

function [param_map, cleanup] = makeBaseParamMap()
test_dir = fileparts(mfilename('fullpath'));
pair_path = fullfile(test_dir, ['norm_pair_', char(java.util.UUID.randomUUID), '.txt']);
fid = fopen(pair_path, 'w');
fprintf(fid, 'protein\tpeptide\nP1\tPEPTIDE\n');
fclose(fid);
cleanup = onCleanup(@() deleteIfExists(pair_path));

param_map = containers.Map('KeyType', 'char', 'ValueType', 'any');
param_map(CPTMdecoderWorkflowParamKeys.PARAM_SPEC_DIR_PATH) = 'data/demo/spectra';
param_map(CPTMdecoderWorkflowParamKeys.PARAM_MS1_TOLERANCE_VALUE) = '10';
param_map(CPTMdecoderWorkflowParamKeys.PARAM_MS1_TOLERANCE_TYPE) = 'ppm';
param_map(CPTMdecoderWorkflowParamKeys.PARAM_ALPHA) = '0.01';
param_map(CPTMdecoderWorkflowParamKeys.PARAM_RESULT_FILTER_THRESHOLD) = '0.1';
param_map(CPTMdecoderWorkflowParamKeys.PARAM_OUTPUT_DIR_PATH) = 'output/demo';
param_map(CPTMdecoderWorkflowParamKeys.PARAM_NORM_PROTEIN_PEPTIDE_PAIR_PATH) = pair_path;
param_map(CPTMdecoderWorkflowParamKeys.PARAM_FILTERED_RES_FILE_PATH) = 'data/demo/filtered_result.txt';
end

function verifyLoggedErrorContains(testCase, func_handle, expected_tag)
did_throw = false;
caught_error = [];
try
    func_handle();
catch me
    did_throw = true;
    caught_error = me;
end
testCase.assertTrue(did_throw, 'Expected function to raise an error.');
testCase.assertEqual(caught_error.identifier, 'CLogger:LoggedError');
testCase.verifyTrue(contains(caught_error.message, expected_tag));
end

function deleteIfExists(path)
if isfile(path)
    delete(path);
end
end
