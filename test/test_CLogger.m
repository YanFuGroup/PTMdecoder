function tests = test_CLogger
% TEST_CLOGGER Unit tests for CLogger facade.
% Input:
%   (none)
% Output:
%   tests (matlab.unittest.Test)

tests = functiontests(localfunctions);
end

function setupOnce(testCase)
% SETUPONCE Save external logger configuration
testCase.TestData.originalConfig = CLogger.getConfig();
end

function teardownOnce(testCase)
% TEARDOWNONCE Restore external logger configuration after ALL tests pass
CLogger.resetForTests();
CLogger.configure(testCase.TestData.originalConfig);
end

function testWriteInfoAndWarn(testCase)
% TESTWRITEINFOANDWARN Verify INFO/WARN entries are written to file.

log_file = fullfile(tempdir, ['ptmdecoder_logger_test_', char(java.util.UUID.randomUUID()), '.log']);
cleanup_log = onCleanup(@() deleteIfExists(log_file));

cfg = struct('enabled', true, 'file_level', 'INFO', 'file_path', log_file, 'to_console', false, 'console_level', 'INFO', 'buffer_size', 1);
CLogger.configure(cfg);
CLogger.info('hello %s', 'world');
CLogger.warn('warn-%d', 1);

testCase.verifyTrue(isfile(log_file));
content = fileread(log_file);
testCase.verifyNotEmpty(regexp(content, '\[INFO\] hello world', 'once'));
testCase.verifyNotEmpty(regexp(content, '\[WARN\] warn-1', 'once'));
end

function testLevelFilterWarnOnly(testCase)
% TESTLEVELFILTERWARNONLY Verify INFO is filtered when level is WARN.

log_file = fullfile(tempdir, ['ptmdecoder_logger_test_', char(java.util.UUID.randomUUID()), '.log']);
cleanup_log = onCleanup(@() deleteIfExists(log_file));

cfg = struct('enabled', true, 'file_level', 'WARN', 'file_path', log_file, 'to_console', false, 'console_level', 'INFO', 'buffer_size', 1);
CLogger.configure(cfg);
CLogger.info('should not appear');
CLogger.warn('should appear');

testCase.verifyTrue(isfile(log_file));
content = fileread(log_file);
testCase.verifyEmpty(regexp(content, 'should not appear', 'once'));
testCase.verifyNotEmpty(regexp(content, 'should appear', 'once'));
end

function testErrorLogsThenThrows(testCase)
% TESTERRORLOGSTHENTHROWS Verify error API logs and then raises exception.

log_file = fullfile(tempdir, ['ptmdecoder_logger_test_', char(java.util.UUID.randomUUID()), '.log']);
cleanup_log = onCleanup(@() deleteIfExists(log_file));

cfg = struct('enabled', true, 'file_level', 'INFO', 'file_path', log_file, 'to_console', false, 'console_level', 'INFO', 'buffer_size', 50);
CLogger.configure(cfg);

testCase.verifyError(@() CLogger.error('fatal-%d', 7), 'CLogger:LoggedError');

content = fileread(log_file);
testCase.verifyNotEmpty(regexp(content, '\[ERROR\] fatal-7', 'once'));
end

function testErrorWithMExceptionRethrowsOriginal(testCase)
% TESTERRORWITHMEXCEPTIONRETHROWSORIGINAL Verify original exception is rethrown.

log_file = fullfile(tempdir, ['ptmdecoder_logger_test_', char(java.util.UUID.randomUUID()), '.log']);
cleanup_log = onCleanup(@() deleteIfExists(log_file));

cfg = struct('enabled', true, 'file_level', 'INFO', 'file_path', log_file, 'to_console', false, 'console_level', 'INFO', 'buffer_size', 50);
CLogger.configure(cfg);

original_me = MException('Demo:Failure', 'demo message');
testCase.verifyError(@() CLogger.error(original_me, 'context-%d', 9), 'Demo:Failure');

content = fileread(log_file);
testCase.verifyNotEmpty(regexp(content, '\[ERROR\] context-9', 'once'));
testCase.verifyNotEmpty(regexp(content, '\[ERROR\] Original exception: \[Demo:Failure\] demo message', 'once'));
end

function testDebugGoesToFileWhenFileLevelIsDebug(testCase)
% TESTDEBUGGOESTOFILEWHENFILELEVELISDEBUG Verify DEBUG is persisted when file level is DEBUG.

log_file = fullfile(tempdir, ['ptmdecoder_logger_test_', char(java.util.UUID.randomUUID()), '.log']);
cleanup_log = onCleanup(@() deleteIfExists(log_file));
cfg = struct('enabled', true, 'file_level', 'DEBUG', 'file_path', log_file, 'to_console', false, 'console_level', 'INFO', 'buffer_size', 1);
CLogger.configure(cfg);
CLogger.debug('debug-msg-%d', 11);

content = fileread(log_file);
testCase.verifyNotEmpty(regexp(content, '\[DEBUG\] debug-msg-11', 'once'));
end

function testBufferedFileWriteNeedsFlush(testCase)
% TESTBUFFEREDFILEWRITENEEDSFLUSH Verify logs are persisted after manual flush under batching mode.

log_file = fullfile(tempdir, ['ptmdecoder_logger_test_', char(java.util.UUID.randomUUID()), '.log']);
cleanup_log = onCleanup(@() deleteIfExists(log_file));

cfg = struct('enabled', true, 'file_level', 'INFO', 'file_path', log_file, 'to_console', false, 'console_level', 'INFO', 'buffer_size', 100);
CLogger.configure(cfg);
CLogger.info('buffered-msg');

if isfile(log_file)
    content_before = fileread(log_file);
    testCase.verifyEmpty(regexp(content_before, 'buffered-msg', 'once'));
end

CLogger.flush();

testCase.verifyTrue(isfile(log_file));
content_after = fileread(log_file);
testCase.verifyNotEmpty(regexp(content_after, '\[INFO\] buffered-msg', 'once'));
end

function deleteIfExists(path)
% DELETEIFEXISTS Delete file if it exists.
if isfile(path)
    delete(path);
end
end
