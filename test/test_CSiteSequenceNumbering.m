function tests = test_CSiteSequenceNumbering
% Regression coverage for sequence-backed numbering in all site outputs.
tests = functiontests(localfunctions);
end

function setup(testCase)
workDir = tempname;
mkdir(workDir);
testCase.TestData.workDir = workDir;
testCase.addTeardown(@() rmdir(workDir, 's'));
oldLogger = CLogger.getConfig();
testCase.addTeardown(@() CLogger.configure(oldLogger));
CLogger.configure(struct('file_path', fullfile(workDir, 'test.log'), ...
    'enabled', true, 'to_console', true, 'console_level', 'WARN'));
writeLines(fullfile(workDir, 'proteins.fasta'), ...
    {'>Pm'; 'MAK'; '>Pa'; 'AK'; '>Empty'; ''; '>InternalM'; 'AMAK'; '>Unmapped'; 'AAAAK'});
writeLines(fullfile(workDir, 'sample.mgf'), ...
    {'BEGIN IONS'; 'TITLE=SpecA'; 'PEPMASS=500.25'; 'CHARGE=2+'; 'END IONS'});
writeLines(fullfile(workDir, 'msms.txt'), {'P\tAK'; 'S\tsample.mgf\tSpecA'; 'AK{Acetyl}\t1'});
end

function testResidueNumberingUsesFullProteinSequence(testCase)
% A peptide beginning with M must not imply that its protein begins with M.
contexts = {'Pm,2;', 'Pa,1;', 'InternalM,2;'};
peptides = {'_AK{Acetyl}_', '_AK{Acetyl}_', '_M{Acetyl}AK_'};
legacyFalse = {'HK2ac', 'AK2ac', 'IM2ac'};
legacyTrue = {'HK3ac', 'AK2ac', 'IM2ac'};
matrixFalse = {'Pm,2; K_ac', 'Pa,2; K_ac', 'InternalM,2; M_ac'};
matrixTrue = {'Pm,3; K_ac', 'Pa,2; K_ac', 'InternalM,2; M_ac'};
for countM = [false, true]
    for idx = 1:numel(contexts)
        cfg = makeConfig(testCase, countM, contexts{idx}, peptides{idx});
        if countM
            expectedLegacy = legacyTrue{idx};
            expectedMatrix = matrixTrue{idx};
        else
            expectedLegacy = legacyFalse{idx};
            expectedMatrix = matrixFalse{idx};
        end
        verifyOutputs(testCase, cfg, {expectedLegacy}, {expectedMatrix});
    end
end
end

function testMultiProteinContextsUseIndependentOffsets(testCase)
cfg = makeConfig(testCase, false, 'Pm,2;Pa,1;', '_AK{Acetyl}_');
verifyOutputs(testCase, cfg, {'HK2ac'}, {'Pm,2;Pa,2; K_ac'});
% Legacy still chooses the first mapped hit, and looks it up by original ID.
cfg = makeConfig(testCase, false, 'Unmapped,4;Pa,1;Pm,2;', '_AK{Acetyl}_');
verifyOutputs(testCase, cfg, {'AK2ac'}, {'Unmapped,5;Pa,2;Pm,2; K_ac'});
end

function testTerminalAndInitialResidueConventionsArePreserved(testCase)
for countM = [false, true]
    cfg = makeConfig(testCase, countM, 'Pm,1;', '_{Acetyl}M{Acetyl}AK_{Acetyl}');
    if countM
        expectedLegacy = {'HN-term_ac', 'HM1ac', 'HC-term_ac'};
        expectedMatrix = {'Pm,0; N-term_ac', 'Pm,1; M_ac', 'Pm,4; C-term_ac'};
    else
        expectedLegacy = {'HN-term_ac', 'HM0ac', 'HC-term_ac'};
        expectedMatrix = {'Pm,-1; N-term_ac', 'Pm,0; M_ac', 'Pm,3; C-term_ac'};
    end
    verifyOutputs(testCase, cfg, expectedLegacy, expectedMatrix);
end
cfg = makeConfig(testCase, false, 'Pm,2;Pa,1;', '_{Acetyl}AK_{Acetyl}');
verifyOutputs(testCase, cfg, {'HN-term_ac', 'HC-term_ac'}, ...
    {'Pm,0;Pa,0; N-term_ac', 'Pm,3;Pa,3; C-term_ac'});
end

function testMissingFastaRaisesConfigurationError(testCase)
cfg = makeConfig(testCase, false, 'Pm,2;Pa,1;', '_AK{Acetyl}_');
for pathValue = {'', "", [], '  '}
    cfg.fasta_file_path = pathValue{1};
    for kind = {'legacy', 'matrix', 'trace'}
        testCase.verifyError(@() runOutput(cfg, kind{1}), 'CSitePositionNumbering:MissingFasta');
    end
end
end

function testMissingAndEmptySequencesRaiseErrorsAndCloseInputs(testCase)
for protein = {'Missing', 'Empty'}
    cfg = makeConfig(testCase, false, [protein{1}, ',1;'], '_AK{Acetyl}_');
    if strcmp(protein{1}, 'Missing')
        expectedId = 'CSitePositionNumbering:MissingProteinSequence';
    else
        expectedId = 'CSitePositionNumbering:EmptyProteinSequence';
    end
    for kind = {'legacy', 'matrix', 'trace'}
        openFiles = fopen('all');
        testCase.verifyError(@() runOutput(cfg, kind{1}), expectedId);
        testCase.verifyEqual(fopen('all'), openFiles);
    end
end
end

function testTerminalOnlyReportsAlsoRequireSequenceEvidence(testCase)
for peptide = {'_{Acetyl}AK_', '_AK_{Acetyl}'}
    for protein = {'Missing', 'Empty'}
        cfg = makeConfig(testCase, false, [protein{1}, ',1;'], peptide{1});
        if strcmp(protein{1}, 'Missing')
            expectedId = 'CSitePositionNumbering:MissingProteinSequence';
        else
            expectedId = 'CSitePositionNumbering:EmptyProteinSequence';
        end
        for kind = {'legacy', 'matrix', 'trace'}
            testCase.verifyError(@() runOutput(cfg, kind{1}), expectedId);
        end
    end
end
end

function testMatrixAndTraceValidateEveryReferencedProtein(testCase)
cfg = makeConfig(testCase, false, 'Pm,2;Missing,1;', '_AK{Acetyl}_');
% The site-list output references only its first mapped protein.
testCase.verifyEqual(runOutput(cfg, 'legacy'), {'HK2ac'});
for kind = {'matrix', 'trace'}
    testCase.verifyError(@() runOutput(cfg, kind{1}), 'CSitePositionNumbering:MissingProteinSequence');
end
end

function testCountingKeepsCoordinatesWithoutSequenceEvidence(testCase)
cfg = makeConfig(testCase, true, 'Missing,1;Empty,1;Pm,2;', '_AK{Acetyl}_');
cfg.fasta_file_path = '';
cfg.regular_express = '[';
verifyOutputs(testCase, cfg, {'UK2ac'}, {'Missing,2;Empty,2;Pm,3; K_ac'});
% A configured FASTA with absent/empty proteins is not consulted either.
cfg.fasta_file_path = fullfile(testCase.TestData.workDir, 'proteins.fasta');
verifyOutputs(testCase, cfg, {'UK2ac'}, {'Missing,2;Empty,2;Pm,3; K_ac'});
end

function testSequenceErrorsIdentifyTheProtein(testCase)
numbering = CSitePositionNumbering(false, ...
    fullfile(testCase.TestData.workDir, 'proteins.fasta'), '^>(\S+)');
for protein = {'Missing', 'Empty'}
    try
        numbering.getOffset(protein{1});
        testCase.verifyFail('Expected a protein sequence error.');
    catch cause
        testCase.verifyTrue(contains(cause.message, protein{1}));
        testCase.verifyTrue(contains(cause.message, 'site_position_count_initial_m=false'));
    end
end
end

function testConfiguredUnavailableFastaFailsButCountingSkipsIt(testCase)
cfg = makeConfig(testCase, false, 'Pm,2;', '_AK{Acetyl}_');
cfg.fasta_file_path = fullfile(testCase.TestData.workDir, 'not_present.fasta');
for kind = {'legacy', 'matrix', 'trace'}
    testCase.verifyError(@() runOutput(cfg, kind{1}), 'CSitePositionNumbering:UnavailableFasta');
end
cfg.site_position_count_initial_m = true;
cfg.regular_express = '[';
verifyOutputs(testCase, cfg, {'HK3ac'}, {'Pm,3; K_ac'});
end

function testInvalidOrMissingRegexFails(testCase)
cfg = makeConfig(testCase, false, 'Pm,2;', '_AK{Acetyl}_');
for regexValue = {'[', '', '(', ')', 'a)', 'a**', '(?invalid)', 'a{2,1}', ...
        '\|([^|]+)\|(', '\|([^|]+)\|)', '\|([^|]+)\|**'}
    cfg.regular_express = regexValue{1};
    for kind = {'legacy', 'matrix', 'trace'}
        testCase.verifyError(@() runOutput(cfg, kind{1}), 'CSitePositionNumbering:InvalidRegularExpression');
    end
end
end

function testLiteralPipeCharacterClassesUseReaderDialect(testCase)
cfg = makeConfig(testCase, false, 'Pm,2;Pa,1;', '_AK{Acetyl}_');
writeLines(cfg.fasta_file_path, {'>sp|Pm|ProteinM'; 'M A'; 'K'; '>sp|Pa|ProteinA'; 'AK'});
expressions = {'\|([^|]+)\|', '\|([^\|]+)\|', ...
    '^>sp[|]([^|]+)[|]', '^>sp\|(?<protein>[^|]+)\|', ...
    string('\|([^|]+)\|')};
maps = cell(size(expressions));
for idx = 1:numel(expressions)
    cfg.regular_express = expressions{idx};
    reader = CFastaReader(cfg.fasta_file_path, cfg.regular_express);
    maps{idx} = reader.read(true);
    testCase.verifyEqual(sort(keys(maps{idx})), {'Pa', 'Pm'});
    testCase.verifyEqual(maps{idx}('Pm'), 'MAK');
    testCase.verifyEqual(maps{idx}('Pa'), 'AK');
    numbering = CSitePositionNumbering(false, cfg.fasta_file_path, cfg.regular_express);
    testCase.verifyEqual(numbering.getOffset('Pm'), 1);
    testCase.verifyEqual(numbering.getOffset('Pa'), 0);
    verifyOutputs(testCase, cfg, {'HK2ac'}, {'Pm,2;Pa,2; K_ac'});
    cfg.site_position_count_initial_m = true;
    verifyOutputs(testCase, cfg, {'HK3ac'}, {'Pm,3;Pa,2; K_ac'});
    cfg.site_position_count_initial_m = false;
end
testCase.verifyEqual(values(maps{1}, {'Pa', 'Pm'}), values(maps{2}, {'Pa', 'Pm'}));
end

function testUpstreamBackreferenceRegexRemainsSupported(testCase)
cfg = makeConfig(testCase, false, 'Pm,2;Pa,1;', '_AK{Acetyl}_');
cfg.regular_express = '^>(\w+)\|\1$';
writeLines(cfg.fasta_file_path, {'>Pm|Pm'; 'MAK'; '>Pa|Pa'; 'AK'});
verifyOutputs(testCase, cfg, {'HK2ac'}, {'Pm,2;Pa,2; K_ac'});
end

function testRegexValidationRestoresCallerWarningState(testCase)
cfg = makeConfig(testCase, false, 'Pm,2;', '_AK{Acetyl}_');
previousState = warning;
[previousBacktrace, previousVerbose] = warningModes();
[previousMessage, previousId] = lastwarn;
cleanup = onCleanup(@() restoreWarnings(previousState, previousMessage, previousId, ...
    previousBacktrace, previousVerbose)); %#ok<NASGU>
warning('off', 'all');
warning('on', 'backtrace');
warning('off', 'verbose');
lastwarn('caller diagnostic', 'test:CallerDiagnostic');
expectedState = warning;
[expectedBacktrace, expectedVerbose] = warningModes();
for expression = {'[', '\|([^|]+)\|(', '\|([^|]+)\|', '^>(\w+)\|\1$'}
    cfg.regular_express = expression{1};
    if strcmp(expression{1}, '[') || endsWith(expression{1}, '(')
        testCase.verifyError(@() CSitePositionNumbering(false, cfg.fasta_file_path, ...
            cfg.regular_express), 'CSitePositionNumbering:InvalidRegularExpression');
    else
        if startsWith(expression{1}, '^>')
            writeLines(cfg.fasta_file_path, {'>Pm|Pm'; 'MAK'});
        else
            writeLines(cfg.fasta_file_path, {'>sp|Pm|ProteinM'; 'MAK'});
        end
        numbering = CSitePositionNumbering(false, cfg.fasta_file_path, cfg.regular_express);
        testCase.verifyEqual(numbering.getOffset('Pm'), 1);
    end
    testCase.verifyEqual(warning, expectedState);
    [backtrace, verbose] = warningModes();
    testCase.verifyEqual(backtrace, expectedBacktrace);
    testCase.verifyEqual(verbose, expectedVerbose);
    [message, id] = lastwarn;
    testCase.verifyEqual(message, 'caller diagnostic');
    testCase.verifyEqual(id, 'test:CallerDiagnostic');
end
end

function [backtrace, verbose] = warningModes()
backtrace = warning('query', 'backtrace');
verbose = warning('query', 'verbose');
end

function restoreWarnings(state, message, id, backtrace, verbose)
warning(state);
warning(backtrace);
warning(verbose);
lastwarn(message, id);
end

function testConflictingFastaEntriesFail(testCase)
cfg = makeConfig(testCase, false, 'Pm,2;', '_AK{Acetyl}_');
% Covers duplicates stored mid-file and at EOF, including differing non-M sequences.
for lines = {{'>Pm'; 'MAK'; '>Pm'; 'AK'}, ...
        {'>Pa'; 'AK'; '>Pa'; 'AAK'; '>Pm'; 'MAK'}}
    writeLines(cfg.fasta_file_path, lines{1});
    for kind = {'legacy', 'matrix', 'trace'}
        testCase.verifyError(@() runOutput(cfg, kind{1}), 'CFastaReader:ConflictingProteinSequence');
    end
end
end

function testIdenticalDuplicatesAndMultilineSequencesAreAccepted(testCase)
cfg = makeConfig(testCase, false, 'Pm,2;Pa,1;', '_AK{Acetyl}_');
writeLines(cfg.fasta_file_path, {'>Pm'; 'M A'; 'K'; '>Pa'; 'AK'; '>Pm'; 'MAK'});
verifyOutputs(testCase, cfg, {'HK2ac'}, {'Pm,2;Pa,2; K_ac'});
end

function testParameterMapsSupplyFastaToAllThreeOutputs(testCase)
cfg = makeConfig(testCase, false, 'Pm,2;Pa,1;', '_AK{Acetyl}_');
params = containers.Map('KeyType', 'char', 'ValueType', 'any');
params('pep_level_file_path') = cfg.input_path;
params('output_intere_path') = cfg.output_intere_path;
params('output_unintere_path') = cfg.output_unintere_path;
params('output_site_dataset_matrix_path') = cfg.output_site_dataset_matrix_path;
params('spec_dir_path') = cfg.spec_dir_path;
params('msms_res_path') = cfg.msms_res_path;
params('output_trace_peptide_msms_path') = cfg.output_trace_peptide_msms_path;
params('output_trace_site_peptide_path') = cfg.output_trace_site_peptide_path;
params('protein_name_abbr_num') = '2';
params('protein_name_abbr_1') = 'Pm>H';
params('protein_name_abbr_2') = 'Pa>A';
params('mod_name_abbr_num') = '1';
params('mod_name_abbr_1') = 'Acetyl>ac';
params('ignore_strings_site_level') = '""';
params('site_position_count_initial_m') = '0';
params('fasta_file_path') = cfg.fasta_file_path;
params('regular_express') = cfg.regular_express;
legacy = CSiteLevelSummary(CSiteLevelPipelineConfig.fromParamMap(params));
legacy = legacy.site_level_summary();
testCase.verifyEqual(keys(legacy.m_result_output_index), {'HK2ac'});
matrix = CSiteLevelDatasetSummary(CSiteLevelDatasetPipelineConfig.fromParamMap(params));
matrix = matrix.site_level_dataset_summary();
testCase.verifyEqual(matrix.m_site_names, {'Pm,2;Pa,2; K_ac'});
trace = CTraceabilityReportService(CTraceabilityReportServiceConfig.fromParamMap(params));
cleanup = onCleanup(@() delete(trace)); %#ok<NASGU>
trace.run();
testCase.verifyEqual(readTraceKeys(cfg.output_trace_site_peptide_path), {'Pm,2;Pa,2; K_ac'});
end

function cfg = makeConfig(testCase, countM, context, peptide)
workDir = testCase.TestData.workDir;
cfg = struct('input_path', fullfile(workDir, 'peptide.txt'), ...
    'output_intere_path', fullfile(workDir, 'site.txt'), ...
    'output_unintere_path', fullfile(workDir, 'uninterested.txt'), ...
    'output_site_dataset_matrix_path', fullfile(workDir, 'matrix.txt'), ...
    'fasta_file_path', fullfile(workDir, 'proteins.fasta'), ...
    'regular_express', '^>(\S+)', 'site_position_count_initial_m', countM, ...
    'ignore_strings', {{}}, 'spec_dir_path', workDir, ...
    'msms_res_path', fullfile(workDir, 'msms.txt'), ...
    'output_trace_peptide_msms_path', fullfile(workDir, 'peptide_trace.txt'), ...
    'output_trace_site_peptide_path', fullfile(workDir, 'site_trace.txt'));
cfg.pep_level_file_path = cfg.input_path;
cfg.protein_name_abbr = containers.Map({'Pm', 'Pa', 'Missing', 'Empty', 'InternalM'}, {'H', 'A', 'U', 'E', 'I'});
cfg.mod_name_abbr = containers.Map({'Acetyl'}, {'ac'});
cfg.column_idxs = struct('icol_seq', 2, 'icol_charge', 3, 'icol_dataset', 4, ...
    'icol_mass_center', 5, 'icol_low_mz_bound', 6, 'icol_high_mz_bound', 7, 'icol_auc', 8);
writeLines(cfg.input_path, {'protein-header'; 'peptide-header'; '@xic-header'; context; ...
    ['*\t', peptide, '\t+2\tsample.mgf\t500.25\t499.95\t500.55\t100']});
end

function verifyOutputs(testCase, cfg, legacyKeys, matrixKeys)
testCase.verifyEqual(runOutput(cfg, 'legacy'), sort(legacyKeys));
testCase.verifyEqual(runOutput(cfg, 'matrix'), sort(matrixKeys));
testCase.verifyEqual(runOutput(cfg, 'trace'), sort(matrixKeys));
end

function siteKeys = runOutput(cfg, kind)
switch kind
    case 'legacy'
        obj = CSiteLevelSummary(CSiteLevelPipelineConfig(cfg));
        obj = obj.site_level_summary();
        obj.write_file();
        siteKeys = sort(keys(obj.m_result_output_index));
    case 'matrix'
        obj = CSiteLevelDatasetSummary(CSiteLevelDatasetPipelineConfig(cfg));
        obj = obj.site_level_dataset_summary();
        obj.write_file();
        siteKeys = obj.m_site_names;
    case 'trace'
        obj = CTraceabilityReportService(cfg);
        cleanup = onCleanup(@() delete(obj)); %#ok<NASGU>
        obj.run();
        siteKeys = readTraceKeys(cfg.output_trace_site_peptide_path);
end
end

function siteKeys = readTraceKeys(path)
lines = regexp(strtrim(fileread(path)), '\r?\n', 'split');
siteKeys = cell(1, numel(lines) - 1);
for idx = 2:numel(lines)
    cols = strsplit(lines{idx}, sprintf('\t'));
    siteKeys{idx - 1} = cols{2};
end
siteKeys = sort(siteKeys);
end

function writeLines(path, lines)
fid = fopen(path, 'w');
cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
for idx = 1:numel(lines)
    fprintf(fid, '%s\n', strrep(lines{idx}, '\t', sprintf('\t')));
end
end
