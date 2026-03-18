classdef test_CIMPProcessingExecutor < matlab.unittest.TestCase
    % Unit tests for CIMPProcessingExecutor

    methods (Test)
        function testBuildBaseGroupsMatchesLegacyAggregation(testCase)
            ms1_tolerance = struct('value', 10, 'isppm', true);
            executor = CIMPProcessingExecutor([], ms1_tolerance, 1, 0.01, 0);

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
    end
end


function state = append_group(state, group)
group.impRtRanges = [];
if isempty(state)
    state = group;
else
    state(end + 1, 1) = group; %#ok<AGROW>
end
end
