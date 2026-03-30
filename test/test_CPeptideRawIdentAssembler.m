function tests = test_CPeptideRawIdentAssembler
% Validate CPeptideRawIdentAssembler MSMS stability filtering behavior.
% Inputs:
%    (none)
% Outputs:
%    tests (matlab.unittest.Test)
%        Function-based test suite.

tests = functiontests(localfunctions);
end


function testFilterByJaccardOnly(testCase)
% Filter one spectrum by Jaccard only.
% Inputs:
%    testCase (matlab.unittest.TestCase)
%        Test case context.
% Outputs:
%    (none)

spectrum_list = [ ...
    buildSpectrum('raw1.mgf', 'scan.100.2.2', 0.55, {'APEPTIDEA'}, 0.6, 0.95, 0.02), ...
    buildSpectrum('raw1.mgf', 'scan.200.2.2', 0.95, {'CPEPTIDEC'}, 0.8, 0.95, 0.02) ...
    ];

deps = makeDeps(struct( ...
    'enabled', true, ...
    'min_jaccard', 0.80, ...
    'min_support_frequency', 0.80, ...
    'max_abundance_mad', 0.10, ...
    'nan_as_fail', true));

[raw_manager, filter_stats] = CPeptideRawIdentAssembler.buildFromSpectrumList(spectrum_list, deps);

testCase.verifyEqual(filter_stats.total_spectra, 2);
testCase.verifyEqual(filter_stats.kept_spectra, 1);
testCase.verifyEqual(filter_stats.dropped_spectra_jaccard, 1);
testCase.verifyEqual(filter_stats.total_imp_candidates, 2);
testCase.verifyEqual(filter_stats.kept_imp_candidates, 1);
testCase.verifyEqual(filter_stats.dropped_imp_candidates, 1);

[raw_names, raw_stores] = raw_manager.getEntries();
testCase.verifyEqual(numel(raw_names), 1);
testCase.verifyEqual(numel(raw_stores), 1);

store = raw_stores{1};
testCase.verifyEqual(store.length, 1);
testCase.verifyEqual(numel(store.impNames), 1);
testCase.verifyEqual(store.impNames{1}, 'CPEPTIDEC');
testCase.verifyEqual(store.ratioMatrix(1, 1), 0.8, 'AbsTol', 1e-12);
end


function testFilterBySupportOnly(testCase)
% Filter one IMP by support frequency only.
% Inputs:
%    testCase (matlab.unittest.TestCase)
%        Test case context.
% Outputs:
%    (none)

spectrum_list = buildSpectrum('raw1.mgf', 'scan.210.2.2', 0.95, ...
    {'APEPTIDEA', 'CPEPTIDEC'}, [0.7, 0.3], [0.85, 0.60], [0.05, 0.05]);

deps = makeDeps(struct( ...
    'enabled', true, ...
    'min_jaccard', 0.80, ...
    'min_support_frequency', 0.80, ...
    'max_abundance_mad', 0.10, ...
    'nan_as_fail', true));

[raw_manager, filter_stats] = CPeptideRawIdentAssembler.buildFromSpectrumList(spectrum_list, deps);

testCase.verifyEqual(filter_stats.total_spectra, 1);
testCase.verifyEqual(filter_stats.kept_spectra, 1);
testCase.verifyEqual(filter_stats.dropped_spectra_jaccard, 0);
testCase.verifyEqual(filter_stats.total_imp_candidates, 2);
testCase.verifyEqual(filter_stats.kept_imp_candidates, 1);
testCase.verifyEqual(filter_stats.dropped_imp_candidates, 1);

[raw_names, raw_stores] = raw_manager.getEntries();
testCase.verifyEqual(numel(raw_names), 1);
store = raw_stores{1};
testCase.verifyEqual(store.length, 1);
testCase.verifyEqual(numel(store.impNames), 1);
testCase.verifyEqual(store.impNames{1}, 'APEPTIDEA');
testCase.verifyEqual(store.ratioMatrix(1, 1), 0.7, 'AbsTol', 1e-12);
end


function testFilterByMadOnly(testCase)
% Filter one IMP by abundance MAD only.
% Inputs:
%    testCase (matlab.unittest.TestCase)
%        Test case context.
% Outputs:
%    (none)

spectrum_list = buildSpectrum('raw1.mgf', 'scan.220.2.2', 0.95, ...
    {'APEPTIDEA', 'CPEPTIDEC'}, [0.65, 0.35], [0.90, 0.90], [0.05, 0.20]);

deps = makeDeps(struct( ...
    'enabled', true, ...
    'min_jaccard', 0.80, ...
    'min_support_frequency', 0.80, ...
    'max_abundance_mad', 0.10, ...
    'nan_as_fail', true));

[raw_manager, filter_stats] = CPeptideRawIdentAssembler.buildFromSpectrumList(spectrum_list, deps);

testCase.verifyEqual(filter_stats.total_spectra, 1);
testCase.verifyEqual(filter_stats.kept_spectra, 1);
testCase.verifyEqual(filter_stats.dropped_spectra_jaccard, 0);
testCase.verifyEqual(filter_stats.total_imp_candidates, 2);
testCase.verifyEqual(filter_stats.kept_imp_candidates, 1);
testCase.verifyEqual(filter_stats.dropped_imp_candidates, 1);

[raw_names, raw_stores] = raw_manager.getEntries();
testCase.verifyEqual(numel(raw_names), 1);
store = raw_stores{1};
testCase.verifyEqual(store.length, 1);
testCase.verifyEqual(numel(store.impNames), 1);
testCase.verifyEqual(store.impNames{1}, 'APEPTIDEA');
testCase.verifyEqual(store.ratioMatrix(1, 1), 0.65, 'AbsTol', 1e-12);
end


function testBuildFromSpectrumListSupportsNanPassMode(testCase)
% Keep entries with NaN metrics when nan_as_fail is disabled.
% Inputs:
%    testCase (matlab.unittest.TestCase)
%        Test case context.
% Outputs:
%    (none)

spectrum_list = buildSpectrum('raw2.mgf', 'scan.300.2.2', NaN, {'APEPTIDEA'}, 1.0, NaN, NaN);

deps = makeDeps(struct( ...
    'enabled', true, ...
    'min_jaccard', 0.80, ...
    'min_support_frequency', 0.80, ...
    'max_abundance_mad', 0.10, ...
    'nan_as_fail', false));

[raw_manager, filter_stats] = CPeptideRawIdentAssembler.buildFromSpectrumList(spectrum_list, deps);

testCase.verifyEqual(filter_stats.kept_spectra, 1);
testCase.verifyEqual(filter_stats.dropped_spectra_jaccard, 0);
testCase.verifyEqual(filter_stats.kept_imp_candidates, 1);

[raw_names, raw_stores] = raw_manager.getEntries();
testCase.verifyEqual(numel(raw_names), 1);
store = raw_stores{1};
testCase.verifyEqual(store.length, 1);
testCase.verifyEqual(store.impNames{1}, 'APEPTIDEA');
end


function testCreateSpectrumListDepsBuildsExpectedFields(testCase)
% Build deps via factory and verify required fields exist.
% Inputs:
%    testCase (matlab.unittest.TestCase)
%        Test case context.
% Outputs:
%    (none)

deps = CPeptideRawIdentAssembler.createSpectrumListDeps( ...
    [], [], [], struct('value', 20, 'isppm', true), cell(0, 3), cell(0, 3), struct('enabled', false));

testCase.verifyTrue(isfield(deps, 'cMgfDatasetIO'));
testCase.verifyTrue(isfield(deps, 'cMs12DatasetIO'));
testCase.verifyTrue(isfield(deps, 'cMsFileMapper'));
testCase.verifyTrue(isfield(deps, 'ms1_tolerance'));
testCase.verifyTrue(isfield(deps, 'fixedModNameMass'));
testCase.verifyTrue(isfield(deps, 'variableModNameMass'));
testCase.verifyTrue(isfield(deps, 'msmsStabilityFilter'));
end


function spectrum = buildSpectrum(dataset_name, spectrum_name, jaccard, imp_strs, imp_abuns, support, mad)
% Build one spectrum entry used by assembler tests.
% Inputs:
%    dataset_name (1 x 1 char/string)
%        Dataset name.
%    spectrum_name (1 x 1 char/string)
%        Spectrum title.
%    jaccard (1 x 1 double)
%        Jaccard stability.
%    imp_strs (1 x N cell)
%        IMP names.
%    imp_abuns (1 x N double)
%        IMP abundances.
%    support (1 x N double)
%        Support frequencies.
%    mad (1 x N double)
%        Abundance MAD values.
% Outputs:
%    spectrum (1 x 1 struct)
%        Spectrum entry.

spectrum = struct( ...
    'dataset_name', dataset_name, ...
    'spectrum_name', spectrum_name, ...
    'jaccard_stability', jaccard, ...
    'peptidoform_list_str', {imp_strs}, ...
    'peptidoform_list_abun', imp_abuns, ...
    'peptidoform_list_support_freq', support, ...
    'peptidoform_list_abundance_mad', mad);
end


function deps = makeDeps(filter_cfg)
% Build assembler dependencies for unit tests.
% Inputs:
%    filter_cfg (1 x 1 struct)
%        Stability filter config.
% Outputs:
%    deps (1 x 1 struct)
%        Dependency struct for buildFromSpectrumList.

deps = struct( ...
    'cMgfDatasetIO', MockAssemblerMgfDatasetIO(), ...
    'cMs12DatasetIO', MockAssemblerMs12DatasetIO(), ...
    'cMsFileMapper', MockAssemblerMsFileMapper(), ...
    'ms1_tolerance', struct('value', 20, 'isppm', true), ...
    'fixedModNameMass', {cell(0, 3)}, ...
    'variableModNameMass', {cell(0, 3)}, ...
    'msmsStabilityFilter', filter_cfg);
end
