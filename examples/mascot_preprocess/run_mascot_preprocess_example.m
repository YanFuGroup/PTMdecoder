function exampleResult = run_mascot_preprocess_example()
% RUN_MASCOT_PREPROCESS_EXAMPLE Run the fictional Mascot preprocess example.
% Input:
%   (none)
% Output:
%   exampleResult (struct)
%       Example-only policy, counts, FDR diagnostics, generated text, and
%       cleanup confirmation. This teaching result is not a stable API.

exampleDir = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(exampleDir));
preprocessDir = fullfile(repoRoot, 'FDR_control_generate_pep_spec_list');
fixturePath = fullfile(exampleDir, 'fixtures', 'minimal_mascot.dat');

originalPath = path;
pathCleanup = onCleanup(@() path(originalPath));
addpath(repoRoot);
addpath(preprocessDir);

if isempty(which('robustfit'))
    error( ...
        'PTMdecoder:MascotPreprocessExample:MissingStatisticsToolbox', ...
        ['This example requires robustfit from Statistics and Machine ' ...
         'Learning Toolbox.']);
end

if ~isfile(fixturePath)
    error( ...
        'PTMdecoder:MascotPreprocessExample:MissingFixture', ...
        'The fictional Mascot fixture is unavailable: %s', fixturePath);
end

outputDir = tempname;
[created, message] = mkdir(outputDir);
if ~created
    error( ...
        'PTMdecoder:MascotPreprocessExample:CreateOutputFailed', ...
        'Failed to create temporary output directory "%s": %s', ...
        outputDir, message);
end
outputCleanup = onCleanup(@() removeDirectoryIfPresent(outputDir));

% All policy values below are fictional and example-only. They are not
% recommended scientific defaults.
tagType = 'Protein';
decoyTag = 'DECOY_';
groupTag = 'EXAMPLE_GROUP';
fdrThreshold = 0.5;
selectedFdrMode = 'GF';

result = ReadDatResult(fixturePath);
[DecoyType, GroupType, ~, scores, numrst, I] = JudgeGroup( ...
    result, tagType, decoyTag, groupTag);
[FDR, Iid, threshold, finalFDR, transferCoefficients, lambdaCoefficients] = ...
    ComputeFDR( ...
        DecoyType, GroupType, scores, numrst, I, fdrThreshold);

groupedResult = result(I(~GroupType));
filteredResult = result(Iid.(selectedFdrMode));

% The example intentionally applies no project-specific post-FDR selector.
% The pepSpec writer expects equal peptides to be contiguous.
selectedResult = filteredResult;
[~, peptideOrder] = sort({selectedResult.peptide});
selectedResult = selectedResult(peptideOrder);

groupResultPath = fullfile(outputDir, 'group_result_mascot.txt');
filteredResultPath = fullfile(outputDir, 'filtered_result_mascot.txt');
pepSpecPath = fullfile(outputDir, 'pepSpecFile.txt');

write_mascot_result_table(groupedResult, groupResultPath);
write_mascot_result_table(filteredResult, filteredResultPath);
write_peptide_spectra_list_file(selectedResult, pepSpecPath);

outputText = struct( ...
    'groupResultMascot', fileread(groupResultPath), ...
    'filteredResultMascot', fileread(filteredResultPath), ...
    'pepSpec', fileread(pepSpecPath));

printOutput('group_result_mascot.txt', outputText.groupResultMascot);
printOutput('filtered_result_mascot.txt', outputText.filteredResultMascot);
printOutput('pepSpecFile.txt', outputText.pepSpec);

% Remove outputs before returning so the example never publishes or leaves
% behind files. The onCleanup object still covers every exceptional path.
rmdir(outputDir, 's');
clear outputCleanup;
if isfolder(outputDir)
    error( ...
        'PTMdecoder:MascotPreprocessExample:CleanupFailed', ...
        'Temporary output directory was not removed: %s', outputDir);
end

exampleResult = struct();
exampleResult.policy = struct( ...
    'tagType', tagType, ...
    'decoyTag', decoyTag, ...
    'groupTag', groupTag, ...
    'fdrThreshold', fdrThreshold, ...
    'selectedFdrMode', selectedFdrMode);
exampleResult.counts = struct( ...
    'parsed', numel(result), ...
    'inGroup', numel(groupedResult), ...
    'filtered', numel(filteredResult));
exampleResult.fdr = struct( ...
    'values', FDR, ...
    'thresholds', threshold, ...
    'finalValues', finalFDR, ...
    'transferCoefficients', transferCoefficients, ...
    'lambdaCoefficients', lambdaCoefficients);
exampleResult.outputText = outputText;
exampleResult.cleanupConfirmed = true;
end


function printOutput(filename, content)
% PRINTOUTPUT Print one generated file with a visible heading.
fprintf('\n=== %s ===\n', filename);
fprintf('%s', content);
if isempty(content) || content(end) ~= newline
    fprintf('\n');
end
end


function removeDirectoryIfPresent(directory)
% REMOVEDIRECTORYIFPRESENT Remove an example-owned temporary directory.
if isfolder(directory)
    rmdir(directory, 's');
end
end
