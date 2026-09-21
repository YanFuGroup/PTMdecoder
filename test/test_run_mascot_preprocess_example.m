function tests = test_run_mascot_preprocess_example
% TEST_RUN_MASCOT_PREPROCESS_EXAMPLE Test the runnable Mascot example.
tests = functiontests(localfunctions);
end


function testRunnableExample(testCase)
% TESTRUNNABLEEXAMPLE Validate the full fictional preprocess call chain.
testDir = fileparts(mfilename('fullpath'));
repoRoot = fileparts(testDir);
exampleDir = fullfile(repoRoot, 'examples', 'mascot_preprocess');

originalPath = path;
pathCleanup = onCleanup(@() path(originalPath));
addpath(exampleDir);
pathBeforeRun = path;

testCase.assumeNotEmpty( ...
    which('robustfit'), ...
    'Statistics and Machine Learning Toolbox is required for this example.');

[consoleText, exampleResult] = evalc('run_mascot_preprocess_example()');

testCase.verifyEqual(path, pathBeforeRun, ...
    'The example must restore every path entry that it changes.');
testCase.verifyTrue(exampleResult.cleanupConfirmed);
testCase.verifyEqual(exampleResult.counts.parsed, 12);
testCase.verifyEqual(exampleResult.counts.inGroup, 10);
testCase.verifyEqual(exampleResult.counts.filtered, 4);
testCase.verifyEqual(exampleResult.policy.tagType, 'Protein');
testCase.verifyEqual(exampleResult.policy.decoyTag, 'DECOY_');
testCase.verifyEqual(exampleResult.policy.groupTag, 'EXAMPLE_GROUP');
testCase.verifyEqual(exampleResult.policy.fdrThreshold, 0.5);
testCase.verifyEqual(exampleResult.policy.selectedFdrMode, 'GF');

verifyFdrDiagnostics(testCase, exampleResult.fdr);
verifyMascotTables(testCase, exampleResult.outputText);
verifyPepSpec(testCase, exampleResult.outputText.pepSpec);

testCase.verifySubstring(consoleText, 'Reading minimal_mascot.dat ... done.');
testCase.verifySubstring(consoleText, '=== group_result_mascot.txt ===');
testCase.verifySubstring(consoleText, '=== filtered_result_mascot.txt ===');
testCase.verifySubstring(consoleText, '=== pepSpecFile.txt ===');

normalizedConsoleText = normalizeNewlines(consoleText);
testCase.verifySubstring( ...
    normalizedConsoleText, ...
    normalizeNewlines(exampleResult.outputText.groupResultMascot));
testCase.verifySubstring( ...
    normalizedConsoleText, ...
    normalizeNewlines(exampleResult.outputText.filteredResultMascot));
testCase.verifySubstring( ...
    normalizedConsoleText, ...
    normalizeNewlines(exampleResult.outputText.pepSpec));
end


function verifyFdrDiagnostics(testCase, fdr)
% VERIFYFDRDIAGNOSTICS Validate every FDR mode and both robust fits.
modeNames = {'GF', 'SF', 'TF'};
expectedFields = {'GF'; 'SF'; 'TF'};
testCase.verifyEqual(fieldnames(fdr.values), expectedFields);
testCase.verifyEqual(fieldnames(fdr.thresholds), expectedFields);
testCase.verifyEqual(fieldnames(fdr.finalValues), expectedFields);

testCase.verifyNumElements(fdr.values.GF, 12);
testCase.verifyNumElements(fdr.values.SF, 10);
testCase.verifyNumElements(fdr.values.TF, 10);
for idxMode = 1:numel(modeNames)
    modeName = modeNames{idxMode};
    testCase.verifyTrue(all(isfinite(fdr.values.(modeName))));
    testCase.verifyTrue(isfinite(fdr.thresholds.(modeName)));
    testCase.verifyTrue(isfinite(fdr.finalValues.(modeName)));
end

testCase.verifyTrue(all(isfinite(fdr.transferCoefficients)));
testCase.verifyGreaterThan(max(abs(fdr.transferCoefficients)), 0);
testCase.verifyTrue(all(isfinite(fdr.lambdaCoefficients)));
end


function verifyMascotTables(testCase, outputText)
% VERIFYMASCOTTABLES Validate both generated 14-column Mascot tables.
expectedHeader = sprintf([ ...
    'Site\tDatasetName\tScan\tSpectrum\tCharge\t' ...
    'Calc_neutral_pepmass\tprecursor_neutral_mass\tmassdiff\t' ...
    'num_match_ions\tpeptide\tprotein\tmodification\t' ...
    'modificationlocation\tScore']);

groupLines = textLines(outputText.groupResultMascot);
filteredLines = textLines(outputText.filteredResultMascot);
testCase.verifyNumElements(groupLines, 11);
testCase.verifyNumElements(filteredLines, 5);
testCase.verifyEqual(groupLines{1}, expectedHeader);
testCase.verifyEqual(filteredLines{1}, expectedHeader);

groupRows = splitTableRows(testCase, groupLines(2:end));
filteredRows = splitTableRows(testCase, filteredLines(2:end));

expectedGroupProteins = { ...
    'EXAMPLE_GROUP_TARGET_01', ...
    'EXAMPLE_GROUP_TARGET_02', ...
    'DECOY_EXAMPLE_GROUP_01', ...
    'EXAMPLE_GROUP_TARGET_03', ...
    'EXAMPLE_GROUP_TARGET_04', ...
    'DECOY_EXAMPLE_GROUP_02', ...
    'DECOY_EXAMPLE_GROUP_03', ...
    'EXAMPLE_GROUP_TARGET_05', ...
    'DECOY_EXAMPLE_GROUP_04', ...
    'EXAMPLE_GROUP_TARGET_06' ...
    };
testCase.verifyEqual(tableColumn(groupRows, 11), expectedGroupProteins);

expectedFilteredProteins = { ...
    'EXAMPLE_GROUP_TARGET_01', ...
    'EXAMPLE_GROUP_TARGET_02', ...
    'EXAMPLE_GROUP_TARGET_03', ...
    'EXAMPLE_GROUP_TARGET_04' ...
    };
testCase.verifyEqual(tableColumn(filteredRows, 11), expectedFilteredProteins);
testCase.verifyEqual( ...
    tableColumn(filteredRows, 3), {'1001', '1002', '1005', '1006'});
testCase.verifyFalse(any(contains(tableColumn(filteredRows, 11), 'DECOY_')));
testCase.verifyFalse(any(contains(tableColumn(groupRows, 11), 'EXAMPLE_OTHER')));

allRows = [groupRows, filteredRows];
testCase.verifyTrue(all(strcmp(tableColumn(allRows, 2), 'fictional_run.mgf')));
testCase.verifyTrue(all(~cellfun('isempty', tableColumn(allRows, 1))));
end


function verifyPepSpec(testCase, content)
% VERIFYPEPSPEC Validate peptide grouping and spectrum order.
expectedLines = { ...
    'ACDEK', ...
    sprintf('fictional_run.mgf\tfictional_run.1001.1001.2'), ...
    sprintf('fictional_run.mgf\tfictional_run.1005.1005.2'), ...
    'FGHIK', ...
    sprintf('fictional_run.mgf\tfictional_run.1002.1002.2'), ...
    'LMNPQR', ...
    sprintf('fictional_run.mgf\tfictional_run.1006.1006.2') ...
    };
testCase.verifyEqual(textLines(content), expectedLines);
end


function rows = splitTableRows(testCase, lines)
% SPLITTABLEROWS Split table rows and enforce the fixed 14-column format.
rows = cell(1, numel(lines));
for idxLine = 1:numel(lines)
    rows{idxLine} = regexp(lines{idxLine}, '\t', 'split');
    testCase.verifyNumElements(rows{idxLine}, 14);
end
end


function values = tableColumn(rows, columnIndex)
% TABLECOLUMN Return one column from split table rows.
values = cellfun( ...
    @(row) row{columnIndex}, rows, 'UniformOutput', false);
end


function lines = textLines(content)
% TEXTLINES Split text independent of platform newline conventions.
content = normalizeNewlines(content);
lines = regexp(content, '\n', 'split');
if ~isempty(lines) && isempty(lines{end})
    lines(end) = [];
end
end


function content = normalizeNewlines(content)
% NORMALIZENEWLINES Convert platform line endings to newline characters.
content = strrep(content, sprintf('\r\n'), newline);
content = strrep(content, sprintf('\r'), newline);
end
