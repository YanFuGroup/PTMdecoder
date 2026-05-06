function tests = test_CMgfDatasetIO
% TEST_CMGFDATASETIO Verify CMgfDatasetIO MGF helpers.
% Input:
%   (none)
% Output:
%   tests (matlab.unittest.Test)

tests = functiontests(localfunctions);
end


function testGetDatasetChargeMapBuildsAndLoadsCache(testCase)
% Verify TITLE and scan-number charge lookup, including cache reload.
% Input:
%   testCase (matlab.unittest.TestCase)
% Output:
%   (none)

workDir = tempname;
mkdir(workDir);
cleanupWork = onCleanup(@() rmdir(workDir, 's'));

datasetName = 'sample.mgf';
mgfPath = fullfile(workDir, datasetName);
writeLines(mgfPath, {
    'BEGIN IONS'
    'TITLE=RawFile.1234.1234.2.0.dta'
    'PEPMASS=500.25'
    'CHARGE=2+'
    '100.1 10'
    'END IONS'
    'BEGIN IONS'
    'TITLE=unparseable-title'
    'PEPMASS=700.25'
    'CHARGE=3+'
    '200.1 20'
    'END IONS'
});

mgfIo = CMgfDatasetIO(workDir);
cleanupIo = onCleanup(@() delete(mgfIo));

chargeMap = mgfIo.get_dataset_charge_map(datasetName);
testCase.verifyTrue(isKey(chargeMap, 'RawFile.1234.1234.2.0.dta'));
testCase.verifyTrue(isKey(chargeMap, '1234'));
testCase.verifyTrue(isKey(chargeMap, 'unparseable-title'));
testCase.verifyEqual(chargeMap('RawFile.1234.1234.2.0.dta'), 2);
testCase.verifyEqual(chargeMap('1234'), 2);
testCase.verifyEqual(chargeMap('unparseable-title'), 3);

cachePath = fullfile(workDir, 'sample_MGF_charge_map.mat');
testCase.verifyTrue(isfile(cachePath));

chargeMapCached = mgfIo.get_dataset_charge_map(datasetName);
testCase.verifyEqual(chargeMapCached('1234'), 2);
testCase.verifyEqual(chargeMapCached('unparseable-title'), 3);

delete(cleanupIo);
delete(cleanupWork);
end


function writeLines(filePath, lines)
% Write test text file with platform-native newlines.
% Input:
%   filePath (char)
%       output file path
%   lines (cellstr)
%       lines to write
% Output:
%   (none)

fid = fopen(filePath, 'w');
if fid == -1
    error('Failed to create test file: %s', filePath);
end
cleanupFile = onCleanup(@() fclose(fid));

for idx = 1:numel(lines)
    fprintf(fid, '%s\n', lines{idx});
end
end
