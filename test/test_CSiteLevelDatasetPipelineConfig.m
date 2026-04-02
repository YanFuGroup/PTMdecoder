function tests = test_CSiteLevelDatasetPipelineConfig
% Validate site-level dataset pipeline config parsing from parameter map.
% Inputs:
%    (none)
% Outputs:
%    tests (matlab.unittest.Test)
%        Function-based test suite.

tests = functiontests(localfunctions);
end


function testFromParamMapParsesRequiredPathsWithoutOutputDir(testCase)
% Verify required input and output paths are parsed directly from params.
% Inputs:
%    testCase (matlab.unittest.TestCase)
%        Test case context.
% Outputs:
%    (none)

param_map = makeBaseParamMap();

cfg = CSiteLevelDatasetPipelineConfig.fromParamMap(param_map);

testCase.verifyEqual(cfg.input_path, 'manual_input.txt');
testCase.verifyEqual(cfg.output_site_dataset_matrix_path, 'manual_output.txt');
end


function testFromParamMapMissingPepLevelPathThrows(testCase)
% Verify missing peptide-level path is rejected.
% Inputs:
%    testCase (matlab.unittest.TestCase)
%        Test case context.
% Outputs:
%    (none)

param_map = makeBaseParamMap();
remove(param_map, CPTMdecoderWorkflowParamKeys.PARAM_PEP_LEVEL_FILE_PATH);

f = @() CSiteLevelDatasetPipelineConfig.fromParamMap(param_map);
testCase.verifyError(f, 'CSiteLevelDatasetPipelineConfig:MissingRequiredParam');
end


function testFromParamMapMissingOutputPathThrows(testCase)
% Verify missing output matrix path is rejected.
% Inputs:
%    testCase (matlab.unittest.TestCase)
%        Test case context.
% Outputs:
%    (none)

param_map = makeBaseParamMap();
remove(param_map, CPTMdecoderWorkflowParamKeys.PARAM_OUTPUT_SITE_DATASET_MATRIX_PATH);

f = @() CSiteLevelDatasetPipelineConfig.fromParamMap(param_map);
testCase.verifyError(f, 'CSiteLevelDatasetPipelineConfig:MissingRequiredParam');
end


function param_map = makeBaseParamMap()
% Build a minimal valid map for site-level dataset config parsing.
% Inputs:
%    (none)
% Outputs:
%    param_map (containers.Map)
%        Valid key-value map for CSiteLevelDatasetPipelineConfig.fromParamMap.

param_map = containers.Map('KeyType', 'char', 'ValueType', 'any');
param_map(CPTMdecoderWorkflowParamKeys.PARAM_PEP_LEVEL_FILE_PATH) = 'manual_input.txt';
param_map(CPTMdecoderWorkflowParamKeys.PARAM_OUTPUT_SITE_DATASET_MATRIX_PATH) = 'manual_output.txt';
param_map(CPTMdecoderWorkflowParamKeys.PARAM_MOD_NAME_ABBR_NUM) = '1';
param_map([CPTMdecoderWorkflowParamKeys.PARAM_PREFIX_MOD_NAME_ABBR, '1']) = 'Acetyl[K]>ac';
param_map(CPTMdecoderWorkflowParamKeys.PARAM_IGNORE_STRINGS_SITE_LEVEL) = '"_"';
end
