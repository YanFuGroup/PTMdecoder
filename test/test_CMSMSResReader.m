function tests = test_CMSMSResReader
% TEST_CMSMSRESREADER Test script for CMSMSResReader IO
% Input:
%   (none)
% Output:
%   tests (matlab.unittest.Test)
    tests = functiontests(localfunctions);
end

function testIO(testCase)
% TESTIO Validate CMSMSResReader read logic
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

reader = CMSMSResReader();
resultObj = reader.read_from_msms_res_file(testFile);

% --- VERIFICATION ---

% 1. Global Structure
% PEPTIDE_C_Empty should be removed because it has no Valid Spectra
% PEPTIDE_B has SpecB1_Empty removed, but SpecB2_Valid remains, so Peptide B remains.
testCase.verifyEqual(length(resultObj.Peptides), 2, ...
    ['Should have 2 peptides (A and B). C should be removed. Got: ' num2str(length(resultObj.Peptides))]);

% 2. Peptide A
p1 = resultObj.Peptides(1);
testCase.verifyTrue(strcmp(p1.peptide_sequence, 'PEPTIDE_A'), 'Peptide A sequence mismatch');
testCase.verifyEqual(length(p1.spectrum_list), 1, 'Peptide A should have 1 spectrum');
testCase.verifyEqual(p1.spectrum_list(1).peptidoform_num, 2, 'SpecA1 count mismatch');
testCase.verifyEqual(p1.spectrum_list(1).peptidoform_list_abun(2), 200, 'Abundance value mismatch');

% 3. Peptide B
p2 = resultObj.Peptides(2);
testCase.verifyTrue(strcmp(p2.peptide_sequence, 'PEPTIDE_B'), 'Peptide B sequence mismatch');
% SpecB1_Empty was empty, so it should be removed. Only SpecB2_Valid remains.
testCase.verifyEqual(length(p2.spectrum_list), 1, 'Peptide B should have 1 spectrum (empty one removed)');
testCase.verifyTrue(strcmp(p2.spectrum_list(1).spectrum_name, 'SpecB2_Valid'), 'Remaining spectrum should be SpecB2_Valid');
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
