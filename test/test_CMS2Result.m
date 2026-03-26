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
testCase.verifyEqual(length(keptSpec.peptidoform_list_support_freq), 60, ...
    'Support-frequency buffer should be trimmed to exact size');
testCase.verifyTrue(all(isnan(keptSpec.peptidoform_list_support_freq)), ...
    'Support-frequency values should default to NaN when not provided');
testCase.verifyEqual(length(keptSpec.peptidoform_list_vif), 60, ...
    'Per-peptidoform VIF buffer should be trimmed to exact size');
testCase.verifyTrue(all(isnan(keptSpec.peptidoform_list_vif)), ...
    'Per-peptidoform VIF values should default to NaN when not provided');
testCase.verifyEqual(length(keptSpec.peptidoform_list_abundance_mad), 60, ...
    'Abundance-MAD buffer should be trimmed to exact size');
testCase.verifyTrue(all(isnan(keptSpec.peptidoform_list_abundance_mad)), ...
    'Abundance-MAD values should default to NaN when not provided');
end

function testMetadataStorage(testCase)
% TESTMETADATASTORAGE Validate spectrum metadata and peptidoform support metadata

res = CMS2Result();
res.addPeptide('PEPTIDE_META');
res.addSpectrum('DatasetM', 'SpecM1', struct( ...
    'jaccard_stability', 0.42, ...
    'vif_all_imp_max', 5.1, ...
    'vif_reported_imp_max', 2.3));
res.addPeptidoform('FormM1', 123.4, struct('support_frequency', 0.7, 'vif', 9.4, 'abundance_mad', 0.05));
res.addPeptidoform('FormM2', 10.0, struct('abundance_mad', 0.2));

res.compress();

spec = res.Peptides(1).spectrum_list(1);
testCase.verifyEqual(spec.jaccard_stability, 0.42, 'AbsTol', 1e-12);
testCase.verifyEqual(spec.vif_all_imp_max, 5.1, 'AbsTol', 1e-12);
testCase.verifyEqual(spec.vif_reported_imp_max, 2.3, 'AbsTol', 1e-12);
testCase.verifyEqual(spec.peptidoform_list_support_freq(1), 0.7, 'AbsTol', 1e-12);
testCase.verifyEqual(spec.peptidoform_list_vif(1), 9.4, 'AbsTol', 1e-12);
testCase.verifyEqual(spec.peptidoform_list_abundance_mad(1), 0.05, 'AbsTol', 1e-12);
testCase.verifyTrue(isnan(spec.peptidoform_list_support_freq(2)), ...
    'Support-frequency should default to NaN when support_frequency field is not provided');
testCase.verifyTrue(isnan(spec.peptidoform_list_vif(2)), ...
    'Per-peptidoform VIF should default to NaN when vif field is not provided');
testCase.verifyEqual(spec.peptidoform_list_abundance_mad(2), 0.2, 'AbsTol', 1e-12);
end


function testAddOrSelectPeptideIndexCache(testCase)
% TESTADDORSELECTPEPTIDEINDEXCACHE Validate map-based peptide lookup behavior

res = CMS2Result();

res.addOrSelectPeptide('PEP_A');
res.addSpectrum('Dataset1', 'SpecA1');
res.addPeptidoform('PEP_A_FORM1', 1.0);

res.addOrSelectPeptide('PEP_B');
res.addSpectrum('Dataset1', 'SpecB1');
res.addPeptidoform('PEP_B_FORM1', 1.0);

% Select existing peptide and append one more spectrum.
res.addOrSelectPeptide('PEP_A');
res.addSpectrum('Dataset1', 'SpecA2');
res.addPeptidoform('PEP_A_FORM2', 1.0);

testCase.verifyEqual(length(res.Peptides), 2, ...
    'Selecting existing peptide must not create duplicate peptide entries.');
testCase.verifyEqual(length(res.Peptides(1).spectrum_list), 2, ...
    'PEP_A should contain two spectra after selecting existing peptide.');
testCase.verifyEqual(length(res.Peptides(2).spectrum_list), 1, ...
    'PEP_B should keep one spectrum.');
testCase.verifyTrue(strcmp(res.Peptides(1).spectrum_list(2).spectrum_name, 'SpecA2'), ...
    'Second spectrum of PEP_A should be SpecA2.');

% After compress, peptide index cache should remain consistent.
res.compress();
res.addOrSelectPeptide('PEP_A');
res.addSpectrum('Dataset1', 'SpecA3');
res.addPeptidoform('PEP_A_FORM3', 1.0);

testCase.verifyEqual(length(res.Peptides), 2, ...
    'compress should not alter peptide count in this case.');
testCase.verifyEqual(length(res.Peptides(1).spectrum_list), 3, ...
    'PEP_A should still be selectable correctly after compress.');
testCase.verifyTrue(strcmp(res.Peptides(1).spectrum_list(3).spectrum_name, 'SpecA3'), ...
    'Third spectrum of PEP_A should be SpecA3 after compress.');
end
