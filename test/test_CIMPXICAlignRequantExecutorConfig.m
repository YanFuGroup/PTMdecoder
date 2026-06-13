function tests = test_CIMPXICAlignRequantExecutorConfig
% Validate CIMPXICAlignRequantExecutorConfig defaults and error contract.
tests = functiontests(localfunctions);
end

function testDefaultsAndExplicitValues(testCase)
base = makeBaseConfig();
cfg_default = CIMPXICAlignRequantExecutorConfig(base);
testCase.verifyEqual(cfg_default.minMSMSnum, 1);
testCase.verifyEqual(cfg_default.resFilterThres, 0);
testCase.verifyEqual(cfg_default.align_options, struct());

base.minMSMSnum = 2;
base.resFilterThres = 0.1;
base.align_options = struct('max_shift', 2);
cfg_explicit = CIMPXICAlignRequantExecutorConfig(base);
testCase.verifyEqual(cfg_explicit.minMSMSnum, 2);
testCase.verifyEqual(cfg_explicit.resFilterThres, 0.1);
testCase.verifyEqual(cfg_explicit.align_options, struct('max_shift', 2));
end

function testInvalidConstructorArgs(testCase)
verifyLoggedErrorContains(testCase, @() CIMPXICAlignRequantExecutorConfig(), ...
    'CIMPXICAlignRequantExecutorConfig:InvalidConstructorArgs');
verifyLoggedErrorContains(testCase, @() CIMPXICAlignRequantExecutorConfig(1), ...
    'CIMPXICAlignRequantExecutorConfig:InvalidConstructorArgs');
end

function testInvalidFieldsUseLoggedErrorContract(testCase)
cases = { ...
    'ms12DatasetIO', [], 'InvalidMs12DatasetIO'; ...
    'ms1_tolerance', struct('value', 10), 'InvalidMs1Tolerance'; ...
    'minMSMSnum', 0, 'InvalidMinMSMSnum'; ...
    'resFilterThres', -1, 'InvalidResFilterThres'; ...
    'aligner', [], 'InvalidAligner'; ...
    'align_strategy', [], 'InvalidAlignStrategy'; ...
    'align_options', 1, 'InvalidAlignOptions'};

for idx = 1:size(cases, 1)
    cfg = makeBaseConfig();
    cfg.(cases{idx, 1}) = cases{idx, 2};
    verifyLoggedErrorContains(testCase, @() CIMPXICAlignRequantExecutorConfig(cfg), ...
        ['CIMPXICAlignRequantExecutorConfig:', cases{idx, 3}]);
end
end

function cfg = makeBaseConfig()
cfg = struct( ...
    'ms12DatasetIO', struct('placeholder', true), ...
    'ms1_tolerance', struct('value', 10, 'isppm', true), ...
    'aligner', struct('placeholder', true), ...
    'align_strategy', 'global');
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
