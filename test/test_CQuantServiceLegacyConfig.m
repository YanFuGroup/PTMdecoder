function tests = test_CQuantServiceLegacyConfig
% Verify quant services preserve the legacy sparse-XIC default.
tests = functiontests(localfunctions);
end

function testMissingAndEmptyThresholdDefaultToFive(testCase)
testCase.verifyEqual(CStructOptionUtils.get(struct(), 'min_xic_nonzero_points', 5), 5);
testCase.verifyEqual(CStructOptionUtils.get( ...
    struct('min_xic_nonzero_points', []), 'min_xic_nonzero_points', 5), 5);
end

function testQuantServicesUseLegacyCompatibleThresholdLookup(testCase)
service_files = { ...
    fullfile('@CPeptideQuantService', 'CPeptideQuantService.m'), ...
    fullfile('@CPeptideAlignRequantService', 'CPeptideAlignRequantService.m'), ...
    fullfile('@CNormalizationQuantService', 'CNormalizationQuantService.m')};

for idx = 1:numel(service_files)
    content = fileread(service_files{idx});
    testCase.verifyTrue(contains(content, ...
        'min_xic_nonzero_points = CStructOptionUtils.get('), service_files{idx});
    testCase.verifyFalse(contains(content, ...
        '''minXicNonzeroPoints'', cfg.min_xic_nonzero_points'), service_files{idx});
    testCase.verifyFalse(contains(content, ...
        '''minXicNonzeroPoints'', msms_cfg.min_xic_nonzero_points'), service_files{idx});
end
end

function testNormalizationQuantUsesConfiguredMinMsmsNum(testCase)
service_file = fullfile('@CNormalizationQuantService', 'CNormalizationQuantService.m');
content = fileread(service_file);

testCase.verifyTrue(contains(content, ...
    'min_msms_num = CStructOptionUtils.get('));
testCase.verifyTrue(contains(content, ...
    'msms_cfg, ''min_MSMS_num'', 1)'));
testCase.verifyTrue(contains(content, ...
    '''minMSMSnum'', min_msms_num'));
testCase.verifyFalse(contains(content, ...
    '''minMSMSnum'', 1'));
end
