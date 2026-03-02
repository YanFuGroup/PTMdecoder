function tests = test_CModificationRegistry
% TEST_CMODIFICATIONREGISTRY Unit tests for centralized modification helper.
% Input:
%   (none)
% Output:
%   tests (matlab.unittest.Test)

tests = functiontests(localfunctions);
end

function testLoadMapFromIniAndCompatibilityWrapper(testCase)
% TESTLOADMAPFROMINIANDCOMPATIBILITYWRAPPER Validate map loading and wrapper compatibility.

ini_path = createTempModifyIni();
cleanup_file = onCleanup(@() deleteIfExists(ini_path)); %#ok<NASGU>

map1 = CModificationRegistry.loadMap(ini_path);
map2 = readModifyInfo(ini_path);

testCase.verifyEqual(map1('Carbamidomethyl[C]'), 57.021464);
testCase.verifyEqual(map1('Oxidation[M]'), 15.994915);
testCase.verifyEqual(map1('null'), 0);
testCase.verifyEqual(map1('NULL'), 0);
testCase.verifyEqual(map1('Null'), 0);
testCase.verifyEqual(map1('unknown'), 0);
testCase.verifyEqual(map2('Carbamidomethyl[C]'), map1('Carbamidomethyl[C]'));
testCase.verifyEqual(map2('Oxidation[M]'), map1('Oxidation[M]'));
end

function testBuildModNameMassAndFromConfig(testCase)
% TESTBUILDMODNAMEMASSANDFROMCONFIG Validate parsing behavior and fromConfig wiring.

ini_path = createTempModifyIni();
cleanup_file = onCleanup(@() deleteIfExists(ini_path)); %#ok<NASGU>

cfg = struct( ...
    'mod_file_path', ini_path, ...
    'fixed_mod', 'Carbamidomethyl[C];', ...
    'variable_mod', 'Oxidation[M];');

[fixedModNameMass, variableModNameMass, mapModification] = CModificationRegistry.fromConfig(cfg);

testCase.verifyEqual(size(fixedModNameMass), [1, 3]);
testCase.verifyEqual(size(variableModNameMass), [1, 3]);
testCase.verifyEqual(fixedModNameMass{1,1}, 'Carbamidomethyl');
testCase.verifyEqual(fixedModNameMass{1,2}, 'C');
testCase.verifyEqual(fixedModNameMass{1,3}, mapModification('Carbamidomethyl[C]'));
testCase.verifyEqual(variableModNameMass{1,1}, 'Oxidation');
testCase.verifyEqual(variableModNameMass{1,2}, 'M');
testCase.verifyEqual(variableModNameMass{1,3}, mapModification('Oxidation[M]'));
end

function ini_path = createTempModifyIni()
% CREATETEMPMODIFYINI Create a minimal temporary modify.ini fixture.

ini_path = [tempname, '.ini'];
fid = fopen(ini_path, 'w');
if fid <= 0
    error('test_CModificationRegistry:CreateIniFailed', ...
        'Failed to create temporary ini file.');
end
cleanup_fid = onCleanup(@() fclose(fid)); %#ok<NASGU>

fprintf(fid, '[modify]\n');
fprintf(fid, 'total=2\n');
fprintf(fid, 'name1=Carbamidomethyl[C]\n');
fprintf(fid, 'Carbamidomethyl[C]=C NORMAL 57.021464 57.051300 0\n');
fprintf(fid, 'name2=Oxidation[M]\n');
fprintf(fid, 'Oxidation[M]=M NORMAL 15.994915 16.000000 0\n');
end

function deleteIfExists(file_path)
% DELETEIFEXISTS Delete file if it exists.

if exist(file_path, 'file')
    delete(file_path);
end
end
