function tests = test_CIMPProcessingExecutorConfig
% Validate CIMPProcessingExecutorConfig defaults and threshold validation.
tests = functiontests(localfunctions);
end

function testDefaultAndExplicitMinXicNonzeroPoints(testCase)
base = makeBaseConfig();
cfg_default = CIMPProcessingExecutorConfig(base);
testCase.verifyEqual(cfg_default.minXicNonzeroPoints, 5);

base.minXicNonzeroPoints = 3;
cfg_explicit = CIMPProcessingExecutorConfig(base);
testCase.verifyEqual(cfg_explicit.minXicNonzeroPoints, 3);
end

function testInvalidMinXicNonzeroPoints(testCase)
base = makeBaseConfig();
base.minXicNonzeroPoints = 0;
testCase.verifyError(@() CIMPProcessingExecutorConfig(base), ...
    'CIMPProcessingExecutorConfig:InvalidMinXicNonzeroPoints');

base.minXicNonzeroPoints = 1.5;
testCase.verifyError(@() CIMPProcessingExecutorConfig(base), ...
    'CIMPProcessingExecutorConfig:InvalidMinXicNonzeroPoints');
end

function cfg = makeBaseConfig()
cfg = struct( ...
    'ms12DatasetIO', struct('placeholder', true), ...
    'ms1_tolerance', struct('value', 10, 'isppm', true), ...
    'minMSMSnum', 1, ...
    'alpha', 0.01, ...
    'resFilterThres', 0.1);
end
