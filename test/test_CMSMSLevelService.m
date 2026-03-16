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
    'case_penalty_intens', 'intens_sum', ...
    'grid_penalty_intens', 'intens_sum', ...
    'case_OLS_intens_weight', 'none', ...
    'n_resamples', 3, ...
    'random_seed', 123, ...
    'relative_threshold', 0.01);

svc = CMSMSLevelService(cfg);
svc.run();

outPath = fullfile(workDir, 'report_msms.txt');
res = CMS2ResultIO.read(outPath);

allJaccard = [];
allSupport = [];
for idxPep = 1:numel(res.Peptides)
    specs = res.Peptides(idxPep).spectrum_list;
    for idxSpec = 1:numel(specs)
        allJaccard(end + 1, 1) = specs(idxSpec).jaccard_stability; %#ok<AGROW>
        nImp = specs(idxSpec).peptidoform_num;
        if nImp > 0
            allSupport = [allSupport; specs(idxSpec).peptidoform_list_support_freq(1:nImp)']; %#ok<AGROW>
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

clear svc res;
delete(cleanupPath);
clear classes;
end


function testRunWithoutStabilityOptionsKeepsNaNAndWarns(testCase)
% Verify run() keeps NaN diagnostics and emits warn when stability options are missing
% Input:
%   testCase (matlab.unittest.TestCase)
% Output:
%   (none)

projectRoot = fileparts(fileparts(mfilename('fullpath')));
mockRoot = fullfile(projectRoot, 'test', 'mocks', 'service_orchestration');
cleanupPath = setupServiceOrchestrationMocks(mockRoot); %#ok<NASGU>

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
cleanupLogger = onCleanup(@() CLogger.resetForTests());

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
for idxPep = 1:numel(res.Peptides)
    specs = res.Peptides(idxPep).spectrum_list;
    for idxSpec = 1:numel(specs)
        allJaccard(end + 1, 1) = specs(idxSpec).jaccard_stability; %#ok<AGROW>
        nImp = specs(idxSpec).peptidoform_num;
        if nImp > 0
            allSupport = [allSupport; specs(idxSpec).peptidoform_list_support_freq(1:nImp)']; %#ok<AGROW>
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

logText = fileread(logPath);
testCase.verifyNotEmpty(strfind(logText, 'Missing stability options'));

clear svc res;
delete(cleanupLogger);
delete(cleanupPath);
clear classes;
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
clear classes;
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
