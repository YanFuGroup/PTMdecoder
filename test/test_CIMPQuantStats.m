function tests = test_CIMPQuantStats
% Validate CIMPQuantStats flush overwrite behavior.
% Inputs:
%    (none)
% Outputs:
%    tests (matlab.unittest.Test)
%        Function-based test suite.

tests = functiontests(localfunctions);
end


function testFlushOverwritesExistingFile(testCase)
% Verify flush overwrites existing stats file instead of appending.
% Inputs:
%    testCase (matlab.unittest.TestCase)
%        Test case context.
% Outputs:
%    (none)

stats_path = fullfile(tempdir, ['rt_stats_overwrite_', char(java.util.UUID.randomUUID), '.mat']);
cleanup = onCleanup(@() deleteIfExists(stats_path));

CIMPQuantStats.rt_sorted_stats('init');
CIMPQuantStats.rt_sorted_stats('record', 10);
CIMPQuantStats.rt_sorted_stats('record', 20);
CIMPQuantStats.rt_sorted_stats('flush', stats_path);

first_data = load(stats_path, 'rt_sorted_counts');
testCase.verifyEqual(first_data.rt_sorted_counts, [10; 20]);

CIMPQuantStats.rt_sorted_stats('init');
CIMPQuantStats.rt_sorted_stats('record', 30);
CIMPQuantStats.rt_sorted_stats('flush', stats_path);

second_data = load(stats_path, 'rt_sorted_counts');
testCase.verifyEqual(second_data.rt_sorted_counts, 30);
end


function deleteIfExists(file_path)
% Delete file if it exists.
% Inputs:
%    file_path (1 x 1 char/string)
%        target path
% Outputs:
%    (none)

if exist(file_path, 'file')
    delete(file_path);
end
end
