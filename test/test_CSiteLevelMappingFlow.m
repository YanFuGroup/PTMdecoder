function tests = test_CSiteLevelMappingFlow
% Unit tests for site-level protein abbreviation mapping flow.
% Inputs:
%    None.
% Outputs:
%    tests (matlab.unittest.Test)
%        Function-based test suite.
tests = functiontests(localfunctions);
end


function testLoadProteinAbbrMapFromTsvConflictThrows(testCase)
% Verify loader throws when one protein maps to different abbreviations.
tsv_path = createTempPath('.tsv');
cleanup_obj = onCleanup(@() deleteIfExists(tsv_path)); %#ok<NASGU>

content_lines = {
    ['ProteinName', char(9), 'Abbr']
    ['P1', char(9), 'H1']
    ['P1', char(9), 'H2']
};
writeTextFile(tsv_path, content_lines);

f = @() CSiteLevelPipelineConfig.loadProteinAbbrMapFromTsv(tsv_path, 'ProteinName', 'Abbr');
testCase.verifyError(f, 'CLogger:LoggedError');
end


function testSiteLevelSummaryUsesFirstResolvableProteinHit(testCase)
% Verify legacy strategy keeps first resolvable protein hit in one protein line.
input_path = createTempPath('.txt');
intere_path = createTempPath('.txt');
unintere_path = createTempPath('.txt');
cleanup_obj = onCleanup(@() cleanupFiles({input_path, intere_path, unintere_path})); %#ok<NASGU>

% Two proteins are both mappable; the first one must win.
content_lines = {
    'header-1'
    'header-2'
    'P1,10;P2,30;'
    '* AK{Acetyl} x x x x x 100'
};
writeTextFile(input_path, content_lines);

protein_map = containers.Map('KeyType', 'char', 'ValueType', 'char');
protein_map('P1') = 'H1';
protein_map('P2') = 'H2';

mod_map = containers.Map('KeyType', 'char', 'ValueType', 'char');
mod_map('Acetyl') = 'ac';

cfg = CSiteLevelPipelineConfig(struct( ...
    'input_path', input_path, ...
    'output_intere_path', intere_path, ...
    'output_unintere_path', unintere_path, ...
    'protein_abbr_input_mode', 'manual', ...
    'protein_name_abbr', protein_map, ...
    'mod_name_abbr', mod_map, ...
    'ignore_strings', {{}}, ...
    'column_idxs', struct('icol_seq', 2, 'icol_auc', 8)));

summary_obj = CSiteLevelSummary(cfg);
summary_obj = summary_obj.site_level_summary();

key_set = keys(summary_obj.m_result_output_index);

testCase.verifyTrue(any(strcmp(key_set, 'H1K9ac')));
testCase.verifyFalse(any(strcmp(key_set, 'H2K29ac')));
end


function testSiteLevelConfigDefaultCountInitialMFalse(testCase)
% Verify site-level config defaults to not counting the initial M.
cfg = CSiteLevelPipelineConfig(struct( ...
    'input_path', 'input.txt', ...
    'output_intere_path', 'interested.txt', ...
    'output_unintere_path', 'uninterested.txt', ...
    'protein_abbr_input_mode', 'manual', ...
    'protein_name_abbr', containers.Map('KeyType', 'char', 'ValueType', 'char'), ...
    'mod_name_abbr', containers.Map('KeyType', 'char', 'ValueType', 'char')));

testCase.verifyFalse(cfg.site_position_count_initial_m);
end


function testSiteLevelSummaryCanCountInitialM(testCase)
% Verify site-level summary can include initial M in site numbering when enabled.
input_path = createTempPath('.txt');
intere_path = createTempPath('.txt');
unintere_path = createTempPath('.txt');
cleanup_obj = onCleanup(@() cleanupFiles({input_path, intere_path, unintere_path})); %#ok<NASGU>

content_lines = {
    'header-1'
    'header-2'
    'P1,10;'
    '* _AK{Acetyl}_ x x x x x 100'
};
writeTextFile(input_path, content_lines);

protein_map = containers.Map('KeyType', 'char', 'ValueType', 'char');
protein_map('P1') = 'H1';

mod_map = containers.Map('KeyType', 'char', 'ValueType', 'char');
mod_map('Acetyl') = 'ac';

cfg = CSiteLevelPipelineConfig(struct( ...
    'input_path', input_path, ...
    'output_intere_path', intere_path, ...
    'output_unintere_path', unintere_path, ...
    'protein_abbr_input_mode', 'manual', ...
    'protein_name_abbr', protein_map, ...
    'mod_name_abbr', mod_map, ...
    'ignore_strings', {{}}, ...
    'site_position_count_initial_m', true, ...
    'column_idxs', struct('icol_seq', 2, 'icol_auc', 8)));

summary_obj = CSiteLevelSummary(cfg);
summary_obj = summary_obj.site_level_summary();

key_set = keys(summary_obj.m_result_output_index);
testCase.verifyTrue(any(strcmp(key_set, 'H1K11ac')));
end


function testSiteLevelDatasetConfigDefaultCountInitialMTrue(testCase)
% Verify site-level dataset config counts initial M by default.
cfg = CSiteLevelDatasetPipelineConfig(struct( ...
    'input_path', 'input.txt', ...
    'output_site_dataset_matrix_path', 'output.txt', ...
    'mod_name_abbr', containers.Map('KeyType', 'char', 'ValueType', 'char')));

testCase.verifyTrue(cfg.site_position_count_initial_m);
end


function testSiteLevelDatasetSummaryUsesAllProteinsInProteinLine(testCase)
% Verify site-level dataset summary keeps all proteins from one protein line.
input_path = createTempPath('.txt');
output_path = createTempPath('.txt');
cleanup_obj = onCleanup(@() cleanupFiles({input_path, output_path})); %#ok<NASGU>

content_lines = {
    'protein-header'
    'peptide-header'
    'xic-header'
    'P1,10;P2,30;'
    ['*', char(9), '_AK{Acetyl}_', char(9), 'x', char(9), 'DS1', char(9), 'x', char(9), 'x', char(9), 'x', char(9), '100']
};
writeTextFile(input_path, content_lines);

mod_map = containers.Map('KeyType', 'char', 'ValueType', 'char');
mod_map('Acetyl') = 'ac';

cfg = CSiteLevelDatasetPipelineConfig(struct( ...
    'input_path', input_path, ...
    'output_site_dataset_matrix_path', output_path, ...
    'mod_name_abbr', mod_map, ...
    'ignore_strings', {{}}, ...
    'column_idxs', struct('icol_seq', 2, 'icol_dataset', 4, 'icol_auc', 8)));

summary_obj = CSiteLevelDatasetSummary(cfg);
summary_obj = summary_obj.site_level_dataset_summary();

site_keys = keys(summary_obj.m_site_dataset_sum);
testCase.verifyTrue(any(strcmp(site_keys, 'P1,11;P2,31; K_ac')));

dataset_sum_map = summary_obj.m_site_dataset_sum('P1,11;P2,31; K_ac');
testCase.verifyTrue(isKey(dataset_sum_map, 'DS1'));
testCase.verifyEqual(dataset_sum_map('DS1'), 100);
end


function testSiteLevelDatasetSummaryCanExcludeInitialM(testCase)
% Verify site-level dataset summary can exclude initial M when disabled.
input_path = createTempPath('.txt');
output_path = createTempPath('.txt');
cleanup_obj = onCleanup(@() cleanupFiles({input_path, output_path})); %#ok<NASGU>

content_lines = {
    'protein-header'
    'peptide-header'
    'xic-header'
    'P1,10;P2,30;'
    ['*', char(9), '_AK{Acetyl}_', char(9), 'x', char(9), 'DS1', char(9), 'x', char(9), 'x', char(9), 'x', char(9), '100']
};
writeTextFile(input_path, content_lines);

mod_map = containers.Map('KeyType', 'char', 'ValueType', 'char');
mod_map('Acetyl') = 'ac';

cfg = CSiteLevelDatasetPipelineConfig(struct( ...
    'input_path', input_path, ...
    'output_site_dataset_matrix_path', output_path, ...
    'mod_name_abbr', mod_map, ...
    'ignore_strings', {{}}, ...
    'site_position_count_initial_m', false, ...
    'column_idxs', struct('icol_seq', 2, 'icol_dataset', 4, 'icol_auc', 8)));

summary_obj = CSiteLevelDatasetSummary(cfg);
summary_obj = summary_obj.site_level_dataset_summary();

site_keys = keys(summary_obj.m_site_dataset_sum);
testCase.verifyTrue(any(strcmp(site_keys, 'P1,10;P2,30; K_ac')));
end


function writeTextFile(file_path, lines)
% Write text lines into file with trailing newline.
fid = fopen(file_path, 'w');
if fid < 0
    error('test_CSiteLevelMappingFlow:FileOpenFailed', 'Cannot open file: %s', file_path);
end
cleanup_obj = onCleanup(@() fclose(fid));
for idx = 1:numel(lines)
    fprintf(fid, '%s\n', lines{idx});
end
end


function path = createTempPath(ext)
% Build a temp file path with random name.
path = fullfile(tempdir, ['ptmdecoder_site_level_test_', char(java.util.UUID.randomUUID()), ext]);
end


function cleanupFiles(file_list)
% Delete all files in list if they exist.
for idx = 1:numel(file_list)
    deleteIfExists(file_list{idx});
end
end


function deleteIfExists(path)
% Delete one file if present.
if isfile(path)
    delete(path);
end
end
