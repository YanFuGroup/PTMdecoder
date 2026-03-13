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
testCase.verifyTrue(all(isnan(p1.spectrum_list(1).peptidoform_list_support_freq)), ...
    'Legacy peptidoform lines should default support_frequency to NaN');

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
src.addSpectrum('Dataset1', 'SpecA1', struct('jaccard_stability', 0.85));
src.addPeptidoform('FormA1', 100, struct('support_frequency', 0.90));
src.addPeptidoform('FormA2', 200, struct('support_frequency', 0.60));

src.addOrSelectPeptide('PEPTIDE_B');
src.addSpectrum('Dataset2', 'SpecB1', struct('jaccard_stability', NaN));
src.addPeptidoform('FormB1', 300.5, struct('support_frequency', 0.75));

src.compress();

% Round-trip: write then read
CMS2ResultIO.write(src, testFile);
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
        testCase.verifyEqual(dstSpec.peptidoform_list_support_freq, srcSpec.peptidoform_list_support_freq, 'AbsTol', 1e-12);
    end
end
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
fprintf(fid, 'S\tDataset1\tSpecA1\tjaccard=0.125000\tunknown_key=abc\n');
fprintf(fid, 'FormA1\t100\tsupport=0.333333\textra=ignored\n');
fprintf(fid, 'FormA2\t200\n');
fclose(fid);

resultObj = CMS2ResultIO.read(testFile);
spec = resultObj.Peptides(1).spectrum_list(1);

testCase.verifyEqual(spec.jaccard_stability, 0.125, 'AbsTol', 1e-12);
testCase.verifyEqual(spec.peptidoform_list_support_freq(1), 0.333333, 'AbsTol', 1e-12);
testCase.verifyTrue(isnan(spec.peptidoform_list_support_freq(2)), ...
    'Missing support field should default to NaN');
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

testCase.verifyError(funcHandle, 'CLogger:LoggedError');

try
    funcHandle();
catch me
    testCase.verifyTrue(contains(me.message, expectedBusinessTag), ...
        ['Expected error message to contain business tag: ', expectedBusinessTag]);
end
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
