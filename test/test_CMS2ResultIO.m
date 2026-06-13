function tests = test_CMS2ResultIO
% TEST_CMS2RESULTIO Test script for CMS2ResultIO read logic
% Input:
%   (none)
% Outputs:
%   tests (matlab.unittest.Test)

tests = functiontests(localfunctions);
end

function testRead(testCase)
% TESTREAD Validate CMS2ResultIO.read logic
% Input:
%   testCase (matlab.unittest.TestCase)
% Output:
%   (none)

% Create dummy file
testFile = fullfile(pwd, 'test_msms_res_temp.txt');
testCase.addTeardown(@() deleteTestFile(testFile)); % Ensure cleanup even if test fails
fid = fopen(testFile, 'w');
if fid < 0
    error('Could not create temp test file.');
end

% File Content:
% P PEPTIDE_A
% S Dataset1 SpecA1
% FormA1 100
% FormA2 200
% P PEPTIDE_B
% S Dataset1 SpecB1_Empty
% S Dataset1 SpecB2_Valid
% FormB1 300
% P PEPTIDE_C_Empty
% (End of file)

fprintf(fid, 'P\tPEPTIDE_A\n');
fprintf(fid, 'S\tDataset1\tSpecA1\n');
fprintf(fid, 'FormA1\t100\n');
fprintf(fid, 'FormA2\t200\n');
fprintf(fid, 'P\tPEPTIDE_B\n');
fprintf(fid, 'S\tDataset1\tSpecB1_Empty\n'); % Empty spectrum, should be removed
fprintf(fid, 'S\tDataset1\tSpecB2_Valid\n');
fprintf(fid, 'FormB1\t300\n');
fprintf(fid, 'P\tPEPTIDE_C_Empty\n'); % Empty Peptide, should be removed
fclose(fid);

resultObj = CMS2ResultIO.read(testFile);

% --- VERIFICATION ---
% 1. Global Structure
testCase.verifyEqual(length(resultObj.Peptides), 2, ...
    ['Should have 2 peptides (A and B). C should be removed. Got: ' num2str(length(resultObj.Peptides))]);

% 2. Peptide A
p1 = resultObj.Peptides(1);
testCase.verifyTrue(strcmp(p1.peptide_sequence, 'PEPTIDE_A'), 'Peptide A sequence mismatch');
testCase.verifyEqual(length(p1.spectrum_list), 1, 'Peptide A should have 1 spectrum');
testCase.verifyEqual(p1.spectrum_list(1).peptidoform_num, 2, 'SpecA1 count mismatch');
testCase.verifyEqual(p1.spectrum_list(1).peptidoform_list_abun(2), 200, 'Abundance value mismatch');
testCase.verifyTrue(isnan(p1.spectrum_list(1).jaccard_stability), ...
    'Legacy S line should default jaccard_stability to NaN');
testCase.verifyTrue(isnan(p1.spectrum_list(1).vif_all_imp_max), ...
    'Legacy S line should default vif_all_imp_max to NaN');
testCase.verifyTrue(isnan(p1.spectrum_list(1).vif_reported_imp_max), ...
    'Legacy S line should default vif_reported_imp_max to NaN');
testCase.verifyTrue(all(isnan(p1.spectrum_list(1).peptidoform_list_support_freq)), ...
    'Legacy peptidoform lines should default support_frequency to NaN');
testCase.verifyTrue(all(isnan(p1.spectrum_list(1).peptidoform_list_vif)), ...
    'Legacy peptidoform lines should default vif to NaN');
testCase.verifyTrue(all(isnan(p1.spectrum_list(1).peptidoform_list_abundance_mad)), ...
    'Legacy peptidoform lines should default abundance_mad to NaN');

% 3. Peptide B
p2 = resultObj.Peptides(2);
testCase.verifyTrue(strcmp(p2.peptide_sequence, 'PEPTIDE_B'), 'Peptide B sequence mismatch');
% SpecB1_Empty was empty, so it should be removed. Only SpecB2_Valid remains.
testCase.verifyEqual(length(p2.spectrum_list), 1, 'Peptide B should have 1 spectrum (empty one removed)');
testCase.verifyTrue(strcmp(p2.spectrum_list(1).spectrum_name, 'SpecB2_Valid'), 'Remaining spectrum should be SpecB2_Valid');
end

function testWriteReadRoundTrip(testCase)
% TESTWRITEREADROUNDTRIP Validate write->read round-trip consistency for CMS2Result
% Input:
%   testCase (matlab.unittest.TestCase)
% Output:
%   (none)

testFile = fullfile(pwd, 'test_msms_res_roundtrip_temp.txt');
testCase.addTeardown(@() deleteTestFile(testFile));

% Build source result object
src = CMS2Result();
src.addOrSelectPeptide('PEPTIDE_A');
src.addSpectrum('Dataset1', 'SpecA1', struct( ...
    'jaccard_stability', 0.85, ...
    'vif_all_imp_max', 4.2, ...
    'vif_reported_imp_max', 1.8));
src.addPeptidoform('FormA1', 100, struct('support_frequency', 0.90, 'vif', 10.5, 'abundance_mad', 0.01));
src.addPeptidoform('FormA2', 200, struct('support_frequency', 0.60, 'vif', 2.4, 'abundance_mad', 0.08));

src.addOrSelectPeptide('PEPTIDE_B');
src.addSpectrum('Dataset2', 'SpecB1', struct( ...
    'jaccard_stability', NaN, ...
    'vif_all_imp_max', NaN, ...
    'vif_reported_imp_max', NaN));
src.addPeptidoform('FormB1', 300.5, struct('support_frequency', 0.75, 'vif', 3.3, 'abundance_mad', 0.12));

src.compress();

% Round-trip: write then read with VIF fields enabled
CMS2ResultIO.write(src, testFile, true);
dst = CMS2ResultIO.read(testFile);

% Verify hierarchical structure is consistent
testCase.verifyEqual(length(dst.Peptides), length(src.Peptides));

for i = 1:length(src.Peptides)
    srcPep = src.Peptides(i);
    dstPep = dst.Peptides(i);
    testCase.verifyEqual(dstPep.peptide_sequence, srcPep.peptide_sequence);
    testCase.verifyEqual(length(dstPep.spectrum_list), length(srcPep.spectrum_list));

    for j = 1:length(srcPep.spectrum_list)
        srcSpec = srcPep.spectrum_list(j);
        dstSpec = dstPep.spectrum_list(j);
        testCase.verifyEqual(dstSpec.dataset_name, srcSpec.dataset_name);
        testCase.verifyEqual(dstSpec.spectrum_name, srcSpec.spectrum_name);
        testCase.verifyEqual(dstSpec.peptidoform_num, srcSpec.peptidoform_num);
        testCase.verifyEqual(dstSpec.peptidoform_list_str, srcSpec.peptidoform_list_str);
        testCase.verifyEqual(dstSpec.peptidoform_list_abun, srcSpec.peptidoform_list_abun, 'AbsTol', 1e-12);
        testCase.verifyEqual(dstSpec.jaccard_stability, srcSpec.jaccard_stability, 'AbsTol', 1e-12);
        testCase.verifyEqual(dstSpec.vif_all_imp_max, srcSpec.vif_all_imp_max, 'AbsTol', 1e-12);
        testCase.verifyEqual(dstSpec.vif_reported_imp_max, srcSpec.vif_reported_imp_max, 'AbsTol', 1e-12);
        testCase.verifyEqual(dstSpec.peptidoform_list_support_freq, srcSpec.peptidoform_list_support_freq, 'AbsTol', 1e-12);
        testCase.verifyEqual(dstSpec.peptidoform_list_vif, srcSpec.peptidoform_list_vif, 'AbsTol', 1e-12);
        testCase.verifyEqual(dstSpec.peptidoform_list_abundance_mad, srcSpec.peptidoform_list_abundance_mad, 'AbsTol', 1e-12);
    end
end
end

function testReadMultipleFilesAndMerge(testCase)
% TESTREADMULTIPLEFILESANDMERGE Validate object-level merge of multiple MS/MS files.
% Input:
%   testCase (matlab.unittest.TestCase)
% Output:
%   (none)

testFile1 = fullfile(pwd, 'test_msms_res_merge_1_temp.txt');
testFile2 = fullfile(pwd, 'test_msms_res_merge_2_temp.txt');
testCase.addTeardown(@() deleteTestFile(testFile1));
testCase.addTeardown(@() deleteTestFile(testFile2));

src1 = CMS2Result();
src1.addOrSelectPeptide('PEPTIDE_SHARED');
src1.addSpectrum('Dataset1', 'SpecA1', struct( ...
    'jaccard_stability', 0.81, ...
    'vif_all_imp_max', 4.1, ...
    'vif_reported_imp_max', 1.7));
src1.addPeptidoform('FormA1', 100, struct('support_frequency', 0.90, 'vif', 7.5, 'abundance_mad', 0.02));
src1.compress();

src2 = CMS2Result();
src2.addOrSelectPeptide('PEPTIDE_SHARED');
src2.addSpectrum('Dataset2', 'SpecA2', struct( ...
    'jaccard_stability', 0.72, ...
    'vif_all_imp_max', 3.2, ...
    'vif_reported_imp_max', 1.2));
src2.addPeptidoform('FormA2', 200, struct('support_frequency', 0.80, 'vif', 6.5, 'abundance_mad', 0.03));
src2.addOrSelectPeptide('PEPTIDE_UNIQUE');
src2.addSpectrum('Dataset2', 'SpecB1');
src2.addPeptidoform('FormB1', 300);
src2.compress();

CMS2ResultIO.write(src1, testFile1, true);
CMS2ResultIO.write(src2, testFile2, true);

merged = CMS2Result.merge({CMS2ResultIO.read(testFile1), CMS2ResultIO.read(testFile2)});

testCase.verifyEqual(length(merged.Peptides), 2);
testCase.verifyEqual(merged.Peptides(1).peptide_sequence, 'PEPTIDE_SHARED');
testCase.verifyEqual(length(merged.Peptides(1).spectrum_list), 2);
testCase.verifyEqual(merged.Peptides(2).peptide_sequence, 'PEPTIDE_UNIQUE');
testCase.verifyEqual(length(merged.Peptides(2).spectrum_list), 1);

sharedSpec1 = merged.Peptides(1).spectrum_list(1);
sharedSpec2 = merged.Peptides(1).spectrum_list(2);
testCase.verifyEqual(sharedSpec1.dataset_name, 'Dataset1');
testCase.verifyEqual(sharedSpec2.dataset_name, 'Dataset2');
testCase.verifyEqual(sharedSpec1.peptidoform_list_support_freq(1), 0.90, 'AbsTol', 1e-12);
testCase.verifyEqual(sharedSpec2.peptidoform_list_vif(1), 6.5, 'AbsTol', 1e-12);
testCase.verifyEqual(sharedSpec2.peptidoform_list_abundance_mad(1), 0.03, 'AbsTol', 1e-12);
end

function testWriteDefaultsToNoVifFields(testCase)
% TESTWRITEDEFAULTSTONOVIFFIELDS Validate default MSMS output omits VIF fields.

testFile = fullfile(pwd, 'test_msms_res_no_vif_temp.txt');
testCase.addTeardown(@() deleteTestFile(testFile));

src = CMS2Result();
src.addOrSelectPeptide('PEPTIDE_A');
src.addSpectrum('Dataset1', 'SpecA1', struct( ...
    'jaccard_stability', 0.85, ...
    'vif_all_imp_max', 4.2, ...
    'vif_reported_imp_max', 1.8));
src.addPeptidoform('FormA1', 100, struct('support_frequency', 0.90, 'vif', 10.5, 'abundance_mad', 0.01));
src.compress();

CMS2ResultIO.write(src, testFile);

content = fileread(testFile);
lines = strsplit(strtrim(content), newline);
testCase.verifyTrue(any(strcmp(lines, sprintf('S\tDataset1\tSpecA1\tjaccard=0.850000'))));
testCase.verifyTrue(any(strcmp(lines, sprintf('FormA1\t100.000000\tsupport=0.900000\tmad=0.010000'))));
testCase.verifyFalse(contains(content, 'vif_all='));
testCase.verifyFalse(contains(content, 'vif_reported='));
testCase.verifyFalse(contains(content, sprintf('\tvif=')));
end

function testReadNamedFields(testCase)
% TESTREADNAMEDFIELDS Validate named field parsing for jaccard/support

testFile = fullfile(pwd, 'test_msms_res_named_fields_temp.txt');
testCase.addTeardown(@() deleteTestFile(testFile));

fid = fopen(testFile, 'w');
if fid < 0
    error('Could not create temp test file.');
end

fprintf(fid, 'P\tPEPTIDE_A\n');
fprintf(fid, 'S\tDataset1\tSpecA1\tjaccard=0.125000\tvif_all=6.500000\tvif_reported=2.250000\tunknown_key=abc\n');
fprintf(fid, 'FormA1\t100\tsupport=0.333333\tvif=11.200000\tmad=0.050000\textra=ignored\n');
fprintf(fid, 'FormA2\t200\n');
fclose(fid);

resultObj = CMS2ResultIO.read(testFile);
spec = resultObj.Peptides(1).spectrum_list(1);

testCase.verifyEqual(spec.jaccard_stability, 0.125, 'AbsTol', 1e-12);
testCase.verifyEqual(spec.vif_all_imp_max, 6.5, 'AbsTol', 1e-12);
testCase.verifyEqual(spec.vif_reported_imp_max, 2.25, 'AbsTol', 1e-12);
testCase.verifyEqual(spec.peptidoform_list_support_freq(1), 0.333333, 'AbsTol', 1e-12);
testCase.verifyEqual(spec.peptidoform_list_vif(1), 11.2, 'AbsTol', 1e-12);
testCase.verifyEqual(spec.peptidoform_list_abundance_mad(1), 0.05, 'AbsTol', 1e-12);
testCase.verifyTrue(isnan(spec.peptidoform_list_support_freq(2)), ...
    'Missing support field should default to NaN');
testCase.verifyTrue(isnan(spec.peptidoform_list_vif(2)), ...
    'Missing vif field should default to NaN');
testCase.verifyTrue(isnan(spec.peptidoform_list_abundance_mad(2)), ...
    'Missing mad field should default to NaN');
end

function testReadInvalidPeptideLineErrorId(testCase)
% TESTREADINVALIDPEPTIDELINEERRORID Validate invalid peptide line is logged with business tag

testFile = fullfile(pwd, 'test_msms_res_invalid_peptide_temp.txt');
testCase.addTeardown(@() deleteTestFile(testFile));

fid = fopen(testFile, 'w');
if fid < 0
    error('Could not create temp test file.');
end
fprintf(fid, 'P\n');
fclose(fid);

verifyLoggedErrorContains(testCase, @() CMS2ResultIO.read(testFile), ...
    '[CMS2ResultIO:InvalidPeptideLine]');
end

function testReadInvalidSpectrumLineErrorId(testCase)
% TESTREADINVALIDSPECTRUMLINEERRORID Validate invalid spectrum line is logged with business tag

testFile = fullfile(pwd, 'test_msms_res_invalid_spectrum_temp.txt');
testCase.addTeardown(@() deleteTestFile(testFile));

fid = fopen(testFile, 'w');
if fid < 0
    error('Could not create temp test file.');
end
fprintf(fid, 'P\tPEPTIDE_A\n');
fprintf(fid, 'S\tDatasetOnly\n');
fclose(fid);

verifyLoggedErrorContains(testCase, @() CMS2ResultIO.read(testFile), ...
    '[CMS2ResultIO:InvalidSpectrumLine]');
end

function testReadBadLineErrorId(testCase)
% TESTREADBADLINEERRORID Validate malformed peptidoform line is logged with business tag

testFile = fullfile(pwd, 'test_msms_res_bad_line_temp.txt');
testCase.addTeardown(@() deleteTestFile(testFile));

fid = fopen(testFile, 'w');
if fid < 0
    error('Could not create temp test file.');
end
fprintf(fid, 'P\tPEPTIDE_A\n');
fprintf(fid, 'S\tDataset1\tSpec1\n');
fprintf(fid, 'BAD\t123\tsupport=abc\n');
fclose(fid);

verifyLoggedErrorContains(testCase, @() CMS2ResultIO.read(testFile), ...
    '[CMS2ResultIO:InvalidPeptidoformNamedField]');
end

function testReadInvalidMadNamedFieldValue(testCase)
% Validate invalid mad named value is logged with business tag

testFile = fullfile(pwd, 'test_msms_res_invalid_mad_named_field_temp.txt');
testCase.addTeardown(@() deleteTestFile(testFile));

fid = fopen(testFile, 'w');
if fid < 0
    error('Could not create temp test file.');
end
fprintf(fid, 'P\tPEPTIDE_A\n');
fprintf(fid, 'S\tDataset1\tSpec1\n');
fprintf(fid, 'IMP_A\t12.3\tmad=abc\n');
fclose(fid);

verifyLoggedErrorContains(testCase, @() CMS2ResultIO.read(testFile), ...
    '[CMS2ResultIO:InvalidPeptidoformNamedField]');
end

function testReadInvalidSpectrumNamedFieldValue(testCase)
% Validate invalid jaccard named value is logged with business tag

testFile = fullfile(pwd, 'test_msms_res_invalid_spectrum_named_field_temp.txt');
testCase.addTeardown(@() deleteTestFile(testFile));

fid = fopen(testFile, 'w');
if fid < 0
    error('Could not create temp test file.');
end
fprintf(fid, 'P\tPEPTIDE_A\n');
fprintf(fid, 'S\tDataset1\tSpec1\tjaccard=abc\n');
fclose(fid);

verifyLoggedErrorContains(testCase, @() CMS2ResultIO.read(testFile), ...
    '[CMS2ResultIO:InvalidSpectrumNamedField]');
end

function testReadEmptyFileErrorId(testCase)
% TESTREADEMPTYFILEERRORID Validate empty file is logged with business tag

testFile = fullfile(pwd, 'test_msms_res_empty_temp.txt');
testCase.addTeardown(@() deleteTestFile(testFile));

fid = fopen(testFile, 'w');
if fid < 0
    error('Could not create temp test file.');
end
fclose(fid);

verifyLoggedErrorContains(testCase, @() CMS2ResultIO.read(testFile), ...
    '[CMS2ResultIO:EmptyFile]');
end

function testReadOpenFileFailedErrorId(testCase)
% TESTREADOPENFILEFAILEDERRORID Validate non-existing path is logged with business tag

missingFile = fullfile(pwd, 'test_msms_res_missing_temp.txt');
verifyLoggedErrorContains(testCase, @() CMS2ResultIO.read(missingFile), ...
    '[CMS2ResultIO:OpenFileFailed]');
end


function verifyLoggedErrorContains(testCase, funcHandle, expectedBusinessTag)
% VERIFYLOGGEDERRORCONTAINS Verify CLogger error id and business tag in message
% Input:
%   testCase (matlab.unittest.TestCase)
%   funcHandle (1 x 1 function_handle)
%       function expected to raise CLogger.error
%   expectedBusinessTag (1 x N char/string)
%       business tag like [CMS2ResultIO:InvalidSpectrumLine]

did_throw = false;
caught_error = [];
try
    funcHandle();
catch me
    did_throw = true;
    caught_error = me;
end
testCase.assertTrue(did_throw, 'Expected function to raise an error.');
testCase.assertEqual(caught_error.identifier, 'CLogger:LoggedError');
testCase.verifyTrue(contains(caught_error.message, expectedBusinessTag), ...
    ['Expected error message to contain business tag: ', expectedBusinessTag]);
end

function deleteTestFile(testFile)
% DELETETESTFILE Remove temp test file if exists
% Input:
%   testFile (1 x N char/string)
% Output:
%   (none)
if exist(testFile, 'file')
    delete(testFile);
end
end
