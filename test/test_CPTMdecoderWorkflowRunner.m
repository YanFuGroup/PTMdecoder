function tests = test_CPTMdecoderWorkflowRunner
% Validate workflow runner stage logs are persisted to file.
% Inputs:
%    (none)
% Outputs:
%    tests (matlab.unittest.Test)
%        Function-based test suite.

tests = functiontests(localfunctions);
end


function setupOnce(testCase)
% Preserve current logger configuration for test isolation.
% Inputs:
%    testCase (matlab.unittest.TestCase)
%        Test case handle.
% Outputs:
%    (none)

testCase.TestData.originalConfig = CLogger.getConfig();
end


function teardownOnce(testCase)
% Restore logger configuration after all tests in this file.
% Inputs:
%    testCase (matlab.unittest.TestCase)
%        Test case handle.
% Outputs:
%    (none)

CLogger.resetForTests();
CLogger.configure(testCase.TestData.originalConfig);
end


function testRunFlushesBufferedStageLogs(testCase)
% Persist stage start/end logs even when file buffer is large.
% Inputs:
%    testCase (matlab.unittest.TestCase)
%        Test case handle.
% Outputs:
%    (none)

log_file = fullfile(tempdir, ['ptmdecoder_workflow_runner_test_', char(java.util.UUID.randomUUID()), '.log']);
cleanup_log = onCleanup(@() deleteIfExists(log_file)); %#ok<NASGU>

cfg = struct( ...
    'enabled', true, ...
    'file_level', 'INFO', ...
    'file_path', log_file, ...
    'to_console', false, ...
    'console_level', 'ERROR', ...
    'buffer_size', 100);
CLogger.configure(cfg);

stage = struct( ...
    'name', CPTMdecoderWorkflowConfig.STAGE_MERGE_TO_PAIR_LEVEL, ...
    'config', {{}}, ...
    'enabled', true);
workflow_cfg = CPTMdecoderWorkflowConfig(struct('stages', {{stage}}));
runner = CPTMdecoderWorkflowRunner(workflow_cfg);
runner.run();

testCase.verifyTrue(isfile(log_file));
content = fileread(log_file);
testCase.verifyNotEmpty(regexp(content, '\[INFO\] Stage start: merge_to_pair_level', 'once'));
testCase.verifyNotEmpty(regexp(content, '\[INFO\] Stage end: merge_to_pair_level \(success\)', 'once'));
end


function deleteIfExists(path)
% Delete file when it exists.
% Inputs:
%    path (1 x N char)
%        File path.
% Outputs:
%    (none)

if isfile(path)
    delete(path);
end
end
