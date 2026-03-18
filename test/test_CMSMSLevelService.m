function tests = test_CMSMSLevelService
% TEST_CMSMSLEVELSERVICE Verify stage orchestration behavior for CMSMSLevelService
% Input:
%   (none)
% Output:
%   tests (matlab.unittest.Test)

tests = functiontests(localfunctions);
end


function testRunBackfillsNonNaNStability(testCase)
% TESTRUNBACKFILLSNONNANSTABILITY Verify run() writes non-NaN jaccard/support after stage-2 orchestration
% Input:
%   testCase (matlab.unittest.TestCase)
% Output:
%   (none)

projectRoot = fileparts(fileparts(mfilename('fullpath')));
mockRoot = fullfile(projectRoot, 'test', 'mocks', 'service_orchestration');
cleanupPath = setupServiceOrchestrationMocks(mockRoot);

workDir = tempname;
mkdir(workDir);
cleanupWork = onCleanup(@() rmdir(workDir, 's'));

pepSpecPath = fullfile(workDir, 'pep_spec.txt');
writeLines(pepSpecPath, {
    'PEPA'
    'rawA	spec1'
    'rawA	spec2'
});

logPath = fullfile(workDir, 'test_logger.log');
CLogger.resetForTests();
CLogger.configure(struct( ...
    'enabled', true, ...
    'to_console', false, ...
    'file_level', 'DEBUG', ...
    'console_level', 'ERROR', ...
    'file_path', logPath, ...
    'buffer_size', 1));
cleanupLogger = onCleanup(@() configureSharedTestLogger());

cfg = struct( ...
    'spec_dir_path', workDir, ...
    'fasta_file_path', fullfile(workDir, 'dummy.fasta'), ...
    'regular_express', '.*', ...
    'filtered_res_file_path', '', ...
    'pep_spec_file_path', pepSpecPath, ...
    'output_dir_path', workDir, ...
    'model', 1, ...
    'method', 1, ...
    'lambda', 0, ...
    'ms1_tolerance', struct('value', 10, 'isppm', true), ...
    'ms2_tolerance', 0.02, ...
    'alpha', 0.2, ...
    'result_filter_threshold', 0.01, ...
    'ion_types', [1, 2], ...
    'enzyme_name', 'Trypsin', ...
    'enzyme_limits', [1, 0], ...
    'min_peptide_length', 1, ...
    'max_peptide_length', 40, ...
    'max_mod_per_peptide', 5, ...
    'case_penalty_intens', 'intens_sum', ...
    'grid_penalty_intens', 'intens_sum', ...
    'case_OLS_intens_weight', 'none', ...
    'stability_options', struct( ...
        'n_resamples', 3, ...
        'random_seed', 123, ...
        'relative_threshold', 0.01));

svc = CMSMSLevelService(cfg);
svc.run();
CLogger.flush();

outPath = fullfile(workDir, 'report_msms.txt');
res = CMS2ResultIO.read(outPath);

allJaccard = [];
allSupport = [];
allMad = [];
for idxPep = 1:numel(res.Peptides)
    specs = res.Peptides(idxPep).spectrum_list;
    for idxSpec = 1:numel(specs)
        allJaccard(end + 1, 1) = specs(idxSpec).jaccard_stability; %#ok<AGROW>
        nImp = specs(idxSpec).peptidoform_num;
        if nImp > 0
            allSupport = [allSupport; specs(idxSpec).peptidoform_list_support_freq(1:nImp)']; %#ok<AGROW>
            allMad = [allMad; specs(idxSpec).peptidoform_list_abundance_mad(1:nImp)']; %#ok<AGROW>
        end
    end
end

testCase.verifyNotEmpty(allJaccard, ...
    'Expected at least one spectrum entry in report_msms output.');
testCase.verifyTrue(any(~isnan(allJaccard)), ...
    'Expected non-NaN jaccard_stability after stage-2 orchestration.');
testCase.verifyNotEmpty(allSupport, ...
    'Expected at least one peptidoform entry in report_msms output.');
testCase.verifyTrue(any(~isnan(allSupport)), ...
    'Expected non-NaN support_frequency after stage-2 orchestration.');
testCase.verifyNotEmpty(allMad, ...
    'Expected at least one peptidoform abundance_mad entry in report_msms output.');
testCase.verifyTrue(any(~isnan(allMad)), ...
    'Expected non-NaN abundance_mad after stage-2 orchestration.');

logText = fileread(logPath);
testCase.verifyNotEmpty(strfind(logText, 'Stage-2 stability started.'));
testCase.verifyNotEmpty(strfind(logText, 'Stage-2 stability done.'));

clear svc res;
delete(cleanupLogger);
delete(cleanupPath);
end


function testRunWithoutStabilityOptionsKeepsNaNAndWarns(testCase)
% Verify run() keeps NaN diagnostics and emits warn when stability options are missing
% Input:
%   testCase (matlab.unittest.TestCase)
% Output:
%   (none)

projectRoot = fileparts(fileparts(mfilename('fullpath')));
mockRoot = fullfile(projectRoot, 'test', 'mocks', 'service_orchestration');
cleanupPath = setupServiceOrchestrationMocks(mockRoot);

workDir = tempname;
mkdir(workDir);
cleanupWork = onCleanup(@() rmdir(workDir, 's'));

pepSpecPath = fullfile(workDir, 'pep_spec.txt');
writeLines(pepSpecPath, {
    'PEPA'
    'rawA	spec1'
    'rawA	spec2'
});

logPath = fullfile(workDir, 'test_logger.log');
CLogger.resetForTests();
CLogger.configure(struct( ...
    'enabled', true, ...
    'to_console', false, ...
    'file_level', 'DEBUG', ...
    'console_level', 'ERROR', ...
    'file_path', logPath, ...
    'buffer_size', 1));
cleanupLogger = onCleanup(@() configureSharedTestLogger());

cfg = struct( ...
    'spec_dir_path', workDir, ...
    'fasta_file_path', fullfile(workDir, 'dummy.fasta'), ...
    'regular_express', '.*', ...
    'filtered_res_file_path', '', ...
    'pep_spec_file_path', pepSpecPath, ...
    'output_dir_path', workDir, ...
    'model', 1, ...
    'method', 1, ...
    'lambda', 0, ...
    'ms1_tolerance', struct('value', 10, 'isppm', true), ...
    'ms2_tolerance', 0.02, ...
    'alpha', 0.2, ...
    'result_filter_threshold', 0.01, ...
    'ion_types', [1, 2], ...
    'enzyme_name', 'Trypsin', ...
    'enzyme_limits', [1, 0], ...
    'min_peptide_length', 1, ...
    'max_peptide_length', 40, ...
    'max_mod_per_peptide', 5, ...
    'case_penalty_intens', 'intens_sum', ...
    'grid_penalty_intens', 'intens_sum', ...
    'case_OLS_intens_weight', 'none');

svc = CMSMSLevelService(cfg);
svc.run();
CLogger.flush();

outPath = fullfile(workDir, 'report_msms.txt');
res = CMS2ResultIO.read(outPath);

allJaccard = [];
allSupport = [];
allMad = [];
for idxPep = 1:numel(res.Peptides)
    specs = res.Peptides(idxPep).spectrum_list;
    for idxSpec = 1:numel(specs)
        allJaccard(end + 1, 1) = specs(idxSpec).jaccard_stability; %#ok<AGROW>
        nImp = specs(idxSpec).peptidoform_num;
        if nImp > 0
            allSupport = [allSupport; specs(idxSpec).peptidoform_list_support_freq(1:nImp)']; %#ok<AGROW>
            allMad = [allMad; specs(idxSpec).peptidoform_list_abundance_mad(1:nImp)']; %#ok<AGROW>
        end
    end
end

testCase.verifyNotEmpty(allJaccard, ...
    'Expected at least one spectrum entry in report_msms output.');
testCase.verifyTrue(all(isnan(allJaccard)), ...
    'Expected all jaccard_stability values to remain NaN without stability options.');
testCase.verifyNotEmpty(allSupport, ...
    'Expected at least one peptidoform entry in report_msms output.');
testCase.verifyTrue(all(isnan(allSupport)), ...
    'Expected all support_frequency values to remain NaN without stability options.');
testCase.verifyNotEmpty(allMad, ...
    'Expected at least one peptidoform abundance_mad entry in report_msms output.');
testCase.verifyTrue(all(isnan(allMad)), ...
    'Expected all abundance_mad values to remain NaN without stability options.');

logText = fileread(logPath);
testCase.verifyNotEmpty(strfind(logText, 'Stage-2 stability is disabled because cfg.stability_options is not configured'));

clear svc res;
delete(cleanupLogger);
delete(cleanupPath);
end


function testRunParallelStage2MatchesSerial(testCase)
% Verify Stage-2 parallel path produces the same output as serial path.
% Input:
%   testCase (matlab.unittest.TestCase)
% Output:
%   (none)

projectRoot = fileparts(fileparts(mfilename('fullpath')));
mockRoot = fullfile(projectRoot, 'test', 'mocks', 'service_orchestration');
cleanupPath = setupServiceOrchestrationMocks(mockRoot);

workDir = tempname;
mkdir(workDir);
cleanupWork = onCleanup(@() rmdir(workDir, 's')); %#ok<NASGU>

pepSpecPath = fullfile(workDir, 'pep_spec.txt');
writeLines(pepSpecPath, {
    'PEPA'
    'rawA	spec1'
    'rawA	spec2'
});

cfgBase = struct( ...
    'spec_dir_path', workDir, ...
    'fasta_file_path', fullfile(workDir, 'dummy.fasta'), ...
    'regular_express', '.*', ...
    'filtered_res_file_path', '', ...
    'pep_spec_file_path', pepSpecPath, ...
    'output_dir_path', workDir, ...
    'model', 1, ...
    'method', 1, ...
    'lambda', 0, ...
    'ms1_tolerance', struct('value', 10, 'isppm', true), ...
    'ms2_tolerance', 0.02, ...
    'alpha', 0.2, ...
    'result_filter_threshold', 0.01, ...
    'ion_types', [1, 2], ...
    'enzyme_name', 'Trypsin', ...
    'enzyme_limits', [1, 0], ...
    'min_peptide_length', 1, ...
    'max_peptide_length', 40, ...
    'max_mod_per_peptide', 5, ...
    'case_penalty_intens', 'intens_sum', ...
    'grid_penalty_intens', 'intens_sum', ...
    'case_OLS_intens_weight', 'none');

cfgSerial = cfgBase;
cfgSerial.stability_options = struct( ...
    'n_resamples', 3, ...
    'random_seed', 123, ...
    'use_parallel', false, ...
    'relative_threshold', 0.01);

serialOutDir = fullfile(workDir, 'serial');
mkdir(serialOutDir);
cfgSerial.output_dir_path = serialOutDir;
svcSerial = CMSMSLevelService(cfgSerial);
svcSerial.run();
serialRes = CMS2ResultIO.read(fullfile(serialOutDir, 'report_msms.txt'));

cfgParallel = cfgBase;
cfgParallel.stability_options = struct( ...
    'n_resamples', 3, ...
    'random_seed', 123, ...
    'use_parallel', true, ...
    'relative_threshold', 0.01);

parallelOutDir = fullfile(workDir, 'parallel');
mkdir(parallelOutDir);
cfgParallel.output_dir_path = parallelOutDir;
svcParallel = CMSMSLevelService(cfgParallel);
svcParallel.run();
parallelRes = CMS2ResultIO.read(fullfile(parallelOutDir, 'report_msms.txt'));

[serialJaccard, serialSupport, serialMad] = collectStabilityColumns(serialRes);
[parallelJaccard, parallelSupport, parallelMad] = collectStabilityColumns(parallelRes);

testCase.verifyEqual(parallelJaccard, serialJaccard, 'AbsTol', 0);
testCase.verifyEqual(parallelSupport, serialSupport, 'AbsTol', 0);
testCase.verifyEqual(parallelMad, serialMad, 'AbsTol', 0);

clear svcSerial svcParallel serialRes parallelRes;
delete(cleanupPath);
end


function testRunFiltersPeptidesByLength(testCase)
% Verify run() filters peptides outside configured length range.
% Input:
%   testCase (matlab.unittest.TestCase)
% Output:
%   (none)

projectRoot = fileparts(fileparts(mfilename('fullpath')));
mockRoot = fullfile(projectRoot, 'test', 'mocks', 'service_orchestration');
cleanupPath = setupServiceOrchestrationMocks(mockRoot);

workDir = tempname;
mkdir(workDir);
cleanupWork = onCleanup(@() rmdir(workDir, 's')); %#ok<NASGU>

pepSpecPath = fullfile(workDir, 'pep_spec.txt');
writeLines(pepSpecPath, {
    'ABC'
    sprintf('rawA\tspec_short')
    'PEPTIDEA'
    sprintf('rawA\tspec_valid')
    'ABCDEFGHIJKLMNOPQRSTUVWXYABCDEFGHIJKLMNOP'
    sprintf('rawA\tspec_long')
});

cfg = struct( ...
    'spec_dir_path', workDir, ...
    'fasta_file_path', fullfile(workDir, 'dummy.fasta'), ...
    'regular_express', '.*', ...
    'filtered_res_file_path', '', ...
    'pep_spec_file_path', pepSpecPath, ...
    'output_dir_path', workDir, ...
    'model', 1, ...
    'method', 1, ...
    'lambda', 0, ...
    'ms1_tolerance', struct('value', 10, 'isppm', true), ...
    'ms2_tolerance', 0.02, ...
    'alpha', 0.2, ...
    'result_filter_threshold', 0.01, ...
    'ion_types', [1, 2], ...
    'enzyme_name', 'Trypsin', ...
    'enzyme_limits', [1, 0], ...
    'min_peptide_length', 7, ...
    'max_peptide_length', 40, ...
    'max_mod_per_peptide', 5, ...
    'case_penalty_intens', 'intens_sum', ...
    'grid_penalty_intens', 'intens_sum', ...
    'case_OLS_intens_weight', 'none');

svc = CMSMSLevelService(cfg);
svc.run();

outPath = fullfile(workDir, 'report_msms.txt');
res = CMS2ResultIO.read(outPath);

testCase.verifyEqual(numel(res.Peptides), 1, ...
    'Only peptides in configured length range should be processed.');
testCase.verifyEqual(res.Peptides(1).peptide_sequence, 'PEPTIDEA');

clear svc res;
delete(cleanupPath);
end


function writeLines(path, lines)
% WRITELINES Write lines to a UTF-8 text file with newline separators
% Input:
%   path (1 x N char/string)
%       Target file path.
%   lines (K x 1 cell)
%       Text lines to write.
% Output:
%   (none)

fid = fopen(path, 'w');
if fid <= 0
    error('test_CMSMSLevelService:CannotOpenFile', 'Cannot open file: %s', path);
end
cleanupFile = onCleanup(@() fclose(fid));
for idx = 1:numel(lines)
    fprintf(fid, '%s\n', lines{idx});
end
end


function cleanupHandle = setupServiceOrchestrationMocks(mockRoot)
% Activate mock classes for service orchestration tests
% Input:
%   mockRoot (1 x N char/string)
%       Root path that contains @class-folder mocks.
% Output:
%   cleanupHandle (1 x 1 onCleanup)
%       Cleanup handle that restores path and class cache.

addpath(mockRoot, '-begin');
clearClassesAndRestoreSharedLogger();
cleanupHandle = onCleanup(@teardownServiceOrchestrationMocks);
end


function teardownServiceOrchestrationMocks()
% Restore path and class cache after mock-based tests
% Output:
%   (none)

projectRoot = fileparts(fileparts(mfilename('fullpath')));
mockRoot = fullfile(projectRoot, 'test', 'mocks', 'service_orchestration');
if contains(path, [mockRoot, pathsep]) || endsWith(path, mockRoot) || startsWith(path, [mockRoot, pathsep])
    rmpath(mockRoot);
end
end


function clearClassesAndRestoreSharedLogger()
% Clear class cache and reconfigure shared test logger
% Output:
%   (none)

clear CMgfDatasetIO
clear CModificationRegistry
clear CMS2QuantSolver
clear CMS2SpectrumPipeline
clear CMsFileMapper
clear CPepProtService
configureSharedTestLogger();
end

function configureSharedTestLogger()
% CONFIGURESHAREDTESTLOGGER Configure shared test logger under test/output
% Output:
%   (none)

projectRoot = fileparts(fileparts(mfilename('fullpath')));
outputDir = fullfile(projectRoot, 'test', 'output');
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

logPath = fullfile(outputDir, 'ptmdecoder_test.log');
cfg = struct();
cfg.file_path = logPath;
CLogger.resetForTests();
CLogger.configure(cfg);
end


function [allJaccard, allSupport, allMad] = collectStabilityColumns(res)
% Collect stability-related columns from CMS2Result into dense vectors.
% Input:
%   res (CMS2Result)
%       Parsed MSMS report object.
% Output:
%   allJaccard (N x 1 double)
%       Per-spectrum jaccard values.
%   allSupport (M x 1 double)
%       Per-peptidoform support frequency values.
%   allMad (M x 1 double)
%       Per-peptidoform abundance MAD values.

allJaccard = [];
allSupport = [];
allMad = [];
for idxPep = 1:numel(res.Peptides)
    specs = res.Peptides(idxPep).spectrum_list;
    for idxSpec = 1:numel(specs)
        allJaccard(end + 1, 1) = specs(idxSpec).jaccard_stability; %#ok<AGROW>
        nImp = specs(idxSpec).peptidoform_num;
        if nImp > 0
            allSupport = [allSupport; specs(idxSpec).peptidoform_list_support_freq(1:nImp)']; %#ok<AGROW>
            allMad = [allMad; specs(idxSpec).peptidoform_list_abundance_mad(1:nImp)']; %#ok<AGROW>
        end
    end
end
end
