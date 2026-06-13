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
