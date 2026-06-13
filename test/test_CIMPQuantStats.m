function tests = test_CIMPQuantStats
% Validate CIMPQuantStats flush overwrite behavior.
% Inputs:
%    (none)
% Outputs:
%    tests (matlab.unittest.Test)
%        Function-based test suite.

tests = functiontests(localfunctions);
end


function setupOnce(testCase)
testCase.TestData.originalLoggerConfig = CLogger.getConfig();
end


function teardownOnce(testCase)
CLogger.resetForTests();
CLogger.configure(testCase.TestData.originalLoggerConfig);
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


function testQuantGroupStatsRecordsOutcomesAndSparseMetrics(testCase)
CIMPQuantStats.quant_group_stats('init', 3);
CIMPQuantStats.quant_group_stats('record', makeDiagnostics('success', 2, 1));
CIMPQuantStats.quant_group_stats('record', makeDiagnostics('sparse_xic_peaks', 1, 1));
CIMPQuantStats.quant_group_stats('record', makeDiagnostics('no_xic_peak', 0, 0));
CIMPQuantStats.quant_group_stats('record', makeDiagnostics('insufficient_psm_inputs', 0, 0));
CIMPQuantStats.quant_group_stats('record', makeDiagnostics('zero_imp_area', 3, 0));

stats = CIMPQuantStats.quant_group_stats('get', []);
testCase.verifyEqual(fieldnames(stats), { ...
    'min_xic_nonzero_points'; 'total_groups'; 'success'; ...
    'insufficient_psm_inputs'; 'no_xic_peak'; 'sparse_xic_peaks'; ...
    'zero_imp_area'; 'unknown'; 'candidate_peaks'; ...
    'sparse_filtered_groups'; 'filtered_sparse_peaks'});
testCase.verifyEqual(stats.min_xic_nonzero_points, 3);
testCase.verifyEqual(stats.total_groups, 5);
testCase.verifyEqual(stats.success, 1);
testCase.verifyEqual(stats.insufficient_psm_inputs, 1);
testCase.verifyEqual(stats.no_xic_peak, 1);
testCase.verifyEqual(stats.sparse_xic_peaks, 1);
testCase.verifyEqual(stats.zero_imp_area, 1);
testCase.verifyEqual(stats.unknown, 0);
testCase.verifyEqual(stats.candidate_peaks, 6);
testCase.verifyEqual(stats.sparse_filtered_groups, 2);
testCase.verifyEqual(stats.filtered_sparse_peaks, 2);
end


function testQuantGroupStatsWarnsOnceAndIgnoresInvalidSparseMetrics(testCase)
log_path = makeTestLogPath('quant_stats_warn');
cleanup = onCleanup(@() deleteIfExists(log_path)); %#ok<NASGU>
configureTestLogger(log_path);

CIMPQuantStats.quant_group_stats('init', 5);
CIMPQuantStats.quant_group_stats('record', makeDiagnostics('future_reason', 1, 2));
CIMPQuantStats.quant_group_stats('record', makeDiagnostics('future_reason', 1, 2));
CIMPQuantStats.quant_group_stats('record', makeDiagnostics('success', -1, 0));
CIMPQuantStats.quant_group_stats('record', makeDiagnostics('success', 1.5, 0));
CIMPQuantStats.quant_group_stats('record', struct('reason', 'zero_imp_area'));
CIMPQuantStats.quant_group_stats('record', struct('reason', 'zero_imp_area'));
CLogger.flush();

stats = CIMPQuantStats.quant_group_stats('get', []);
testCase.verifyEqual(stats.total_groups, 6);
testCase.verifyEqual(stats.success, 2);
testCase.verifyEqual(stats.zero_imp_area, 2);
testCase.verifyEqual(stats.unknown, 2);
testCase.verifyEqual(stats.candidate_peaks, 0);
testCase.verifyEqual(stats.sparse_filtered_groups, 0);
testCase.verifyEqual(stats.filtered_sparse_peaks, 0);

content = fileread(log_path);
testCase.verifyEqual(count(content, 'CIMPQuantStats:UnknownDiagnosticsReason'), 1);
testCase.verifyEqual(count(content, 'CIMPQuantStats:InvalidSparseMetrics'), 1);
testCase.verifyEqual(count(content, 'CIMPQuantStats:MissingSparseMetrics'), 1);
end


function testLogQuantGroupSummaryWritesThresholdAndAllFields(testCase)
log_path = makeTestLogPath('quant_stats_summary');
cleanup = onCleanup(@() deleteIfExists(log_path)); %#ok<NASGU>
configureTestLogger(log_path);

CIMPQuantStats.quant_group_stats('init', 3);
CIMPQuantStats.quant_group_stats('record', makeDiagnostics('success', 2, 1));
CIMPQuantStats.log_quant_group_summary('Peptide-level');
CLogger.flush();

content = fileread(log_path);
testCase.verifyTrue(contains(content, 'Peptide-level quant group summary:'));
testCase.verifyTrue(contains(content, 'min_xic_nonzero_points=3'));
testCase.verifyTrue(contains(content, 'total_groups=1'));
testCase.verifyTrue(contains(content, 'candidate_peaks=2'));
testCase.verifyTrue(contains(content, 'sparse_filtered_groups=1'));
testCase.verifyTrue(contains(content, 'filtered_sparse_peaks=1'));
end


function diagnostics = makeDiagnostics(reason, candidate_count, filtered_count)
diagnostics = struct( ...
    'reason', reason, ...
    'candidate_peak_count', candidate_count, ...
    'filtered_sparse_peak_count', filtered_count);
end


function configureTestLogger(log_path)
CLogger.resetForTests();
CLogger.configure(struct( ...
    'enabled', true, ...
    'file_level', 'DEBUG', ...
    'file_path', log_path, ...
    'to_console', false, ...
    'console_level', 'INFO', ...
    'buffer_size', 1));
end


function log_path = makeTestLogPath(prefix)
test_dir = fileparts(mfilename('fullpath'));
log_path = fullfile(test_dir, [prefix, '_', char(java.util.UUID.randomUUID), '.log']);
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
