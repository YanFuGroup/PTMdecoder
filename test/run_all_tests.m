function run_all_tests()
% RUN_ALL_TESTS Run all unit tests in the 'test' directory

import matlab.unittest.TestSuite;
import matlab.unittest.TestRunner;
import matlab.unittest.plugins.DiagnosticsRecordingPlugin;

clc;
testDir = fileparts(mfilename('fullpath'));

% Add root code directory to path so helper classes are found
projectDir = fileparts(testDir); 
addpath(projectDir);

fprintf('===================================================\n');
fprintf('Running All Unit Tests\n');
fprintf('===================================================\n\n');

% Create a suite from the folder
suite = TestSuite.fromFolder(testDir);

% Debug: Check if suite is empty
if isempty(suite)
    fprintf('No tests found in folder: %s\n', testDir);
    files = dir(fullfile(testDir, 'test_*.m'));
    if isempty(files)
        fprintf('No files matching test_*.m found.\n');
    else
        fprintf('Found files: %s\n', strjoin({files.name}, ', '));
    end
    return;
end
    
% Create a runner
runner = TestRunner.withTextOutput();

% Run!
result = runner.run(suite);

% Display summary
summary_table = table(result);
disp(summary_table);

% Summarize results
numFailed = sum([result.Failed]);
if numFailed == 0
    fprintf('\nAll tests passed successfully!\n');
else
    fprintf('\nSome tests failed. Number of failed tests: %d\n', numFailed);
end
end