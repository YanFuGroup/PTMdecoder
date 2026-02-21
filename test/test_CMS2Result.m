function tests = test_CMS2Result
% TEST_CMS2RESULT Test script for CMS2Result logic
% Outputs:
%   tests (matlab.unittest.Test)

tests = functiontests(localfunctions);
end

function testLogic(testCase)
% TESTLOGIC Validate CMS2Result add/compress logic
% Input:
%   testCase (matlab.unittest.TestCase)
% Output:
%   (none)

% Use verify functions instead of assert for better reporting
import matlab.unittest.constraints.IsEqualTo;

res = CMS2Result();

% --- Case 1.1: Basic Structure ---
% Add Peptide 1
res.addPeptide('PEPTIDEONE');
% Add Spectrum 1.1
res.addSpectrum('Dataset1', 'Spec1');
res.addPeptidoform('PEPTIDE{Mod}ONE', 1000);
res.addPeptidoform('PEPTIDEONE', 2000);

% --- Case 1.2: Empty Spectrum trimming ---
% Add Peptide 2
res.addPeptide('PEPTIDETWO');
% Add Spectrum 2.1 (Empty, expects to be removed by compress)
res.addSpectrum('Dataset1', 'Spec2FromPep2');

% --- Case 1.3: Buffer Expansion ---
% Add Spectrum 2.2 (Large number of forms to trigger buffer expansion)
res.addSpectrum('Dataset1', 'SpecBufferTest');
% Default buffer is usually 50 in code, let's add 60
for i = 1:60
    res.addPeptidoform(['Form', num2str(i)], i);
end

% Note: Before compress(), length() returns the BUFFERED size, not the logical size.

testCase.verifyEqual(length(res.Peptides), 2, 'Should have 2 peptides initialized');
% Verify buffering behavior (Design Feature): length should be >= logical count
testCase.verifyGreaterThanOrEqual(length(res.Peptides(1).spectrum_list), 1, 'Peptide 1 spectrum buffer size >= 1');
testCase.verifyEqual(res.Peptides(1).spectrum_list(1).peptidoform_num, 2, 'Spec 1 should have 2 peptidoforms');
% Peptide 2 should have 2 logical spectra. Physical length might be 20.
testCase.verifyGreaterThanOrEqual(length(res.Peptides(2).spectrum_list), 2, 'Peptide 2 spectrum buffer size >= 2');

% --- ACTION: Compress ---
res.compress();

% --- VERIFICATION ---
% Verify Peptide 1
testCase.verifyEqual(length(res.Peptides(1).spectrum_list), 1, 'Peptide 1 should still have 1 spectrum');
% Verify Peptide 2
% Spec2FromPep2 was empty, so it should be removed.
% SpecBufferTest has data, should remain.
testCase.verifyEqual(length(res.Peptides(2).spectrum_list), 1, ...
    'Peptide 2 should have 1 spectrum after compress (empty one removed)');

keptSpec = res.Peptides(2).spectrum_list(1);
testCase.verifyTrue(strcmp(keptSpec.spectrum_name, 'SpecBufferTest'), ...
    ['Remaining spectrum should be SpecBufferTest. Actual: ' keptSpec.spectrum_name]);

testCase.verifyEqual(length(keptSpec.peptidoform_list_abun), 60, ...
    ['Buffer should be trimmed to exact size. Expected 60, got ' num2str(length(keptSpec.peptidoform_list_abun))]);
testCase.verifyEqual(keptSpec.peptidoform_list_abun(60), 60, 'Data integrity check failed for last element');
end
