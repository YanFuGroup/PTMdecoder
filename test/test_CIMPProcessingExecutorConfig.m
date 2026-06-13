function tests = test_CIMPProcessingExecutorConfig
% Validate CIMPProcessingExecutorConfig defaults and error contract.
tests = functiontests(localfunctions);
end

function testDefaultsAndExplicitValues(testCase)
base = makeBaseConfig();
cfg_default = CIMPProcessingExecutorConfig(base);
testCase.verifyEqual(cfg_default.minMSMSnum, 1);
testCase.verifyEqual(cfg_default.alpha, 0);
testCase.verifyEqual(cfg_default.resFilterThres, 0);
testCase.verifyEqual(cfg_default.minXicNonzeroPoints, 5);

base.minMSMSnum = 2;
base.alpha = 0.01;
base.resFilterThres = 0.1;
base.minXicNonzeroPoints = 3;
cfg_explicit = CIMPProcessingExecutorConfig(base);
testCase.verifyEqual(cfg_explicit.minMSMSnum, 2);
testCase.verifyEqual(cfg_explicit.alpha, 0.01);
testCase.verifyEqual(cfg_explicit.resFilterThres, 0.1);
testCase.verifyEqual(cfg_explicit.minXicNonzeroPoints, 3);
end

function testInvalidConstructorArgs(testCase)
verifyLoggedErrorContains(testCase, @() CIMPProcessingExecutorConfig(), ...
    'CIMPProcessingExecutorConfig:InvalidConstructorArgs');
verifyLoggedErrorContains(testCase, @() CIMPProcessingExecutorConfig(1), ...
    'CIMPProcessingExecutorConfig:InvalidConstructorArgs');
end

function testInvalidFieldsUseLoggedErrorContract(testCase)
cases = { ...
    'ms12DatasetIO', [], 'InvalidMs12DatasetIO'; ...
    'ms1_tolerance', struct('value', 10), 'InvalidMs1Tolerance'; ...
    'minMSMSnum', 0, 'InvalidMinMSMSnum'; ...
    'alpha', -1, 'InvalidAlpha'; ...
    'resFilterThres', -1, 'InvalidResFilterThres'};

for idx = 1:size(cases, 1)
    cfg = makeBaseConfig();
    cfg.(cases{idx, 1}) = cases{idx, 2};
    verifyLoggedErrorContains(testCase, @() CIMPProcessingExecutorConfig(cfg), ...
        ['CIMPProcessingExecutorConfig:', cases{idx, 3}]);
end
end

function testInvalidMinXicNonzeroPointsUsesLoggedErrorContract(testCase)
invalid_values = {0, -1, 1.5, NaN, Inf};
for idx = 1:numel(invalid_values)
    cfg = makeBaseConfig();
    cfg.minXicNonzeroPoints = invalid_values{idx};
    verifyLoggedErrorContains(testCase, @() CIMPProcessingExecutorConfig(cfg), ...
        'CIMPProcessingExecutorConfig:InvalidMinXicNonzeroPoints');
end
end

function cfg = makeBaseConfig()
cfg = struct( ...
    'ms12DatasetIO', struct('placeholder', true), ...
    'ms1_tolerance', struct('value', 10, 'isppm', true));
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
