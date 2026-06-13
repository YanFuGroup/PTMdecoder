classdef test_CIMPProcessingExecutor < matlab.unittest.TestCase
    % Unit tests for CIMPProcessingExecutor

    methods (Test)
        function testBuildBaseGroupsMatchesLegacyAggregation(testCase)
            ms1_tolerance = struct('value', 10, 'isppm', true);
            executor = makeExecutorWithConfig(struct('placeholder', true), ms1_tolerance);

            rawIdentManager = CIMPRawIdentManager();
            rawStore = rawIdentManager.getOrCreate('rawA.mgf');
            rawStore = rawStore.appendSpecQuant(1.0, 1200, 500.25, 2, ...
                {'PEP{Acetyl[K]}A', 'PEP{Methyl[K]}A'}, [1000.0, 1000.003], [0.7, 0.3]);
            rawStore = rawStore.appendSpecQuant(1.2, 900, 500.25, 2, ...
                {'PEP{Acetyl[K]}A', 'PEP{Methyl[K]}A'}, [1000.0, 1000.003], [0.6, 0.4]);
            rawStore = rawStore.appendSpecQuant(2.0, 800, 600.20, 3, ...
                {'PEP{Trimethyl[K]}A'}, 1200.005, 1.0);
            rawIdentManager.setStore('rawA.mgf', rawStore);

            base_groups = executor.buildBaseGroups(rawIdentManager);

            [raw_names, raw_ident_stores] = rawIdentManager.getEntries();
            aggregator = CIMPGroupAggregator(ms1_tolerance);
            legacy_groups = CIMPGroup.empty(0, 1);
            legacy_groups = aggregator.aggregate(raw_names, raw_ident_stores, [], ...
                @(state, group) append_group(state, group), legacy_groups);

            testCase.verifyEqual(numel(base_groups), numel(legacy_groups));
            for idx = 1:numel(base_groups)
                testCase.verifyEqual(base_groups(idx).rawName, legacy_groups(idx).rawName);
                testCase.verifyEqual(base_groups(idx).selectedCharge, legacy_groups(idx).selectedCharge);
                testCase.verifyEqual(base_groups(idx).impNames, legacy_groups(idx).impNames);
                testCase.verifyEqual(base_groups(idx).impMass, legacy_groups(idx).impMass, 'AbsTol', 1e-10);
                testCase.verifyEqual(base_groups(idx).chargeGroupIdxs, legacy_groups(idx).chargeGroupIdxs);
            end
        end

        function testLegacyConstructorWarnsAndRemainsUsable(testCase)
            ms1_tolerance = struct('value', 10, 'isppm', true);
            executor = [];

            testCase.verifyWarning(@constructLegacyExecutor, ...
                'CIMPProcessingExecutor:DeprecatedLegacyConstructor');
            base_groups = executor.buildBaseGroups([]);

            testCase.verifyEmpty(base_groups);

            function constructLegacyExecutor()
                executor = CIMPProcessingExecutor([], ms1_tolerance, 1, 0.01, 0);
            end
        end

        function testConfigConstructorDoesNotWarn(testCase)
            ms1_tolerance = struct('value', 10, 'isppm', true);
            testCase.verifyWarningFree(@() makeExecutorWithConfig( ...
                struct('placeholder', true), ms1_tolerance));
        end

        function testSparseFilterDebugDistinguishesPartialAndAllFiltered(testCase)
            log_path = makeTestLogPath('executor_sparse_debug');
            cleanup = configureTestLogger(log_path); %#ok<NASGU>
            CIMPQuantStats.quant_group_stats('init', 5);

            partial_executor = makeExecutor(makeDataset([0; 10; 20; 10; zeros(7, 1); 10; 20; 20; 10; 5; zeros(5, 1)]));
            partial_group = makeGroup([0.2; 1.3]);
            partial_executor.quantifyPeptideBlock({'P1', 1}, CIMPRawIdentManager(), partial_group);

            all_executor = makeExecutor(makeDataset([0; 10; 20; 10; zeros(17, 1)]));
            all_group = makeGroup(0.2);
            all_executor.quantifyPeptideBlock({'P1', 1}, CIMPRawIdentManager(), all_group);
            CLogger.flush();

            content = fileread(log_path);
            testCase.verifyTrue(contains(content, 'Partial sparse peak filter'));
            testCase.verifyTrue(contains(content, 'All sparse peaks filtered'));
            testCase.verifyTrue(contains(content, 'candidate_peak_count=2'));
            testCase.verifyTrue(contains(content, 'filtered_sparse_peak_count=1'));
            testCase.verifyTrue(contains(content, 'candidate_nonzero_points=['));
            testCase.verifyFalse(contains(content, 'min_xic_nonzero_points='));
        end

        function testNonSparseFailureDoesNotWriteGroupDebug(testCase)
            log_path = makeTestLogPath('executor_no_sparse_debug');
            cleanup = configureTestLogger(log_path); %#ok<NASGU>
            CIMPQuantStats.quant_group_stats('init', 5);

            executor = makeExecutor(makeDataset(zeros(21, 1)));
            group = makeGroup(0.2);
            executor.quantifyPeptideBlock({'P1', 1}, CIMPRawIdentManager(), group);
            CLogger.flush();

            content = fileread(log_path);
            testCase.verifyFalse(contains(content, '[CIMPProcessingExecutor:onGroupQuant]'));
            stats = CIMPQuantStats.quant_group_stats('get', []);
            testCase.verifyEqual(stats.no_xic_peak, 1);
        end
    end
end


function executor = makeExecutorWithConfig(dataset, ms1_tolerance)
cfg = CIMPProcessingExecutorConfig(struct( ...
    'ms12DatasetIO', dataset, ...
    'ms1_tolerance', ms1_tolerance, ...
    'minMSMSnum', 1, ...
    'alpha', 0.01, ...
    'resFilterThres', 0));
executor = CIMPProcessingExecutor(cfg);
end


function executor = makeExecutor(dataset)
cfg = CIMPProcessingExecutorConfig(struct( ...
    'ms12DatasetIO', dataset, ...
    'ms1_tolerance', struct('value', 10, 'isppm', true), ...
    'minMSMSnum', 1, ...
    'alpha', 0.01, ...
    'resFilterThres', 0, ...
    'minXicNonzeroPoints', 5));
executor = CIMPProcessingExecutor(cfg);
end


function dataset = makeDataset(mono_intensity)
num_scans = numel(mono_intensity);
ms1_index = [(1:num_scans)', (0:0.1:(num_scans - 1) * 0.1)', (2:num_scans + 1)'];
ms1_peaks = [500 * ones(num_scans, 1), mono_intensity];

index_map = containers.Map();
peaks_map = containers.Map();
index_map('raw') = ms1_index;
peaks_map('raw') = ms1_peaks;

file_mapper = struct();
file_mapper.get_ms1_stem = @(name) name;
dataset = struct( ...
    'm_cMsFileMapper', file_mapper, ...
    'm_mapNameMS1Index', index_map, ...
    'm_mapNameMS1Peaks', peaks_map);
end


function group = makeGroup(rts)
num_rows = numel(rts);
group = CIMPGroup( ...
    'raw.mgf', {'IMP_A'}, 1000, ones(num_rows, 1), rts, ...
    100 * ones(num_rows, 1), 2 * ones(num_rows, 1), ...
    499.9, 500.1, 2, (1:num_rows)', []);
end


function cleanup = configureTestLogger(log_path)
original_config = CLogger.getConfig();
CLogger.resetForTests();
CLogger.configure(struct( ...
    'enabled', true, ...
    'file_level', 'DEBUG', ...
    'file_path', log_path, ...
    'to_console', false, ...
    'console_level', 'INFO', ...
    'buffer_size', 1));
cleanup = onCleanup(@() restoreLoggerAndDelete(original_config, log_path));
end


function restoreLoggerAndDelete(original_config, log_path)
CLogger.resetForTests();
CLogger.configure(original_config);
if isfile(log_path)
    delete(log_path);
end
end


function log_path = makeTestLogPath(prefix)
test_dir = fileparts(mfilename('fullpath'));
log_path = fullfile(test_dir, [prefix, '_', char(java.util.UUID.randomUUID), '.log']);
end


function state = append_group(state, group)
group.impRtRanges = [];
if isempty(state)
    state = group;
else
    state(end + 1, 1) = group; %#ok<AGROW>
end
end
