function tests = test_CLoggerCore
% TEST_CLOGGERCORE Unit tests for CLoggerCore behavior.
% Input:
%   (none)
% Output:
%   tests (matlab.unittest.Test)

tests = functiontests(localfunctions);
end

function testConfigureAcceptsLogicalFlags(testCase)
% TESTCONFIGUREACCEPTSLOGICALFLAGS Verify logical config values are accepted.

log_file = fullfile(tempdir, ['ptmdecoder_core_test_', char(java.util.UUID.randomUUID()), '.log']);
cleanup_log = onCleanup(@() deleteIfExists(log_file)); %#ok<NASGU>

core = CLoggerCore();
core.configure(struct( ...
    'enabled', true, ...
    'to_console', false, ...
    'file_level', 'INFO', ...
    'console_level', 'ERROR', ...
    'file_path', log_file, ...
    'buffer_size', 1));

core.log('INFO', 'core-config-ok');
core.flush();

testCase.verifyTrue(isfile(log_file));
content = fileread(log_file);
testCase.verifyNotEmpty(regexp(content, '\[INFO\] core-config-ok', 'once'));
end

function testConfigureRejectsNonLogicalFlags(testCase)
% TESTCONFIGUREREJECTSNONLOGICALFLAGS Verify non-logical config values are rejected.

core = CLoggerCore();
testCase.verifyError(@() core.configure(struct('enabled', '1')), ...
    'CLoggerCore:InvalidLogicalOption');
end

function testBufferingRequiresManualFlush(testCase)
% TESTBUFFERINGREQUIRESMANUALFLUSH Verify buffered entries persist after flush.

log_file = fullfile(tempdir, ['ptmdecoder_core_test_', char(java.util.UUID.randomUUID()), '.log']);
cleanup_log = onCleanup(@() deleteIfExists(log_file)); %#ok<NASGU>

core = CLoggerCore();
core.configure(struct( ...
    'enabled', true, ...
    'to_console', false, ...
    'file_level', 'DEBUG', ...
    'console_level', 'INFO', ...
    'file_path', log_file, ...
    'buffer_size', 100));

core.log('INFO', 'buffered-core-msg');
if isfile(log_file)
    content_before = fileread(log_file);
    testCase.verifyEmpty(regexp(content_before, 'buffered-core-msg', 'once'));
end

core.flush();

testCase.verifyTrue(isfile(log_file));
content_after = fileread(log_file);
testCase.verifyNotEmpty(regexp(content_after, '\[INFO\] buffered-core-msg', 'once'));
end

function testErrorForcesImmediateFlush(testCase)
% TESTERRORFORCESIMMEDIATEFLUSH Verify ERROR entries are flushed immediately.

log_file = fullfile(tempdir, ['ptmdecoder_core_test_', char(java.util.UUID.randomUUID()), '.log']);
cleanup_log = onCleanup(@() deleteIfExists(log_file)); %#ok<NASGU>

core = CLoggerCore();
core.configure(struct( ...
    'enabled', true, ...
    'to_console', false, ...
    'file_level', 'INFO', ...
    'console_level', 'ERROR', ...
    'file_path', log_file, ...
    'buffer_size', 100));

core.log('ERROR', 'core-error-%d', 3);

testCase.verifyTrue(isfile(log_file));
content = fileread(log_file);
testCase.verifyNotEmpty(regexp(content, '\[ERROR\] core-error-3', 'once'));
end

function deleteIfExists(path)
% DELETEIFEXISTS Delete file if it exists.
if isfile(path)
    delete(path);
end
end
