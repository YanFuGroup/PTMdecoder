function tests = test_CPepProtService
% TEST_CPEPPROTSERVICE Unit tests for CPepProtService
% Input:
%   (none)
% Output:
%   tests (matlab.unittest.Test)
tests = functiontests(localfunctions);
end


function setupOnce(testCase)
% SETUPONCE Prepare shared test data
% Input:
%   testCase (matlab.unittest.TestCase)
% Output:
%   (none)
% Create Fasta File
fastaContent = ['>sp|P12345|ProtA' newline ...
                'MKAWVLK' newline ...
                '>sp|P67890|ProtB' newline ...
                'XXXXX' newline ...
                '>sp|P_MULTI_1|Desc1' newline ...
                'MULTI' newline ...
                '>sp|P_MULTI_2|Desc2' newline ...
                'XMULTI'];

testCase.TestData.fastaParams.filename = 'temp_service_test.fasta';
fid = fopen(testCase.TestData.fastaParams.filename, 'w');
fprintf(fid, '%s', fastaContent);
fclose(fid);

% Create Filtered Res File (Tab-delimited)
% Columns: peptide, protein
resContent = sprintf('peptide\tprotein\nMKAWVLK\tP12345\nKAW\tP12345\nWVLK\tP12345\nMULTI\tP_MULTI_1,P_MULTI_2\nXXXXX\tP67890');
                
testCase.TestData.resParams.filename = 'temp_service_res.txt';
fid = fopen(testCase.TestData.resParams.filename, 'w');
fprintf(fid, '%s', resContent);
fclose(fid);

% Initialize Service
% Regex to capture Pxxxxx
regex = '>sp\|([^\|]+)\|';
testCase.TestData.service = CPepProtService(testCase.TestData.fastaParams.filename, ...
    regex, testCase.TestData.resParams.filename);                                        
end


function teardownOnce(testCase)
% TEARDOWNONCE Clean up shared test data
% Input:
%   testCase (matlab.unittest.TestCase)
% Output:
%   (none)
delete(testCase.TestData.fastaParams.filename);
delete(testCase.TestData.resParams.filename);
end


function testGetProteinNamePos_Single(testCase)
% TESTGETPROTEINNAMEPOS_SINGLE Validate single protein hit
% Input:
%   testCase (matlab.unittest.TestCase)
% Output:
%   (none)
srv = testCase.TestData.service;

% Test MKAWVLK in P12345
res = srv.get_protein_name_pos('MKAWVLK');
% Expect {'P12345', 1}
verifyEqual(testCase, size(res, 1), 1);
verifyEqual(testCase, res{1,1}, 'P12345');
verifyEqual(testCase, res{1,2}, 1);

% Test KAW in P12345 (Pos 2)
res = srv.get_protein_name_pos('KAW');
verifyEqual(testCase, size(res, 1), 1);
verifyEqual(testCase, res{1,1}, 'P12345');
verifyEqual(testCase, res{1,2}, 2);
end


function testGetProteinNamePos_Multi(testCase)
% TESTGETPROTEINNAMEPOS_MULTI Validate multiple protein hits
% Input:
%   testCase (matlab.unittest.TestCase)
% Output:
%   (none)
srv = testCase.TestData.service;

% Test MULTI in P_MULTI_1 (Pos 1) and P_MULTI_2 (Pos 2)
res = srv.get_protein_name_pos('MULTI');

verifyEqual(testCase, size(res, 1), 2);

% Order might depend on file read order or hashing. Check content.
names = res(:, 1);

verifyTrue(testCase, any(strcmp(names, 'P_MULTI_1')));
verifyTrue(testCase, any(strcmp(names, 'P_MULTI_2')));

verifyEqual(testCase, res{strcmp(names, 'P_MULTI_1'), 2}, 1);

verifyEqual(testCase, res{strcmp(names, 'P_MULTI_2'), 2}, 2);
end


function testGetWhetherProtNC(testCase)
% TESTGETWHETHERPROTNC Validate N/C terminus flags
% Input:
%   testCase (matlab.unittest.TestCase)
% Output:
%   (none)
srv = testCase.TestData.service;

% MKAWVLK: Pos 1, Len 7, Total 7 -> N=True, C=True
[isN, isC] = srv.getWhetherProtNC('MKAWVLK');
verifyTrue(testCase, isN, 'MKAWVLK should be N-term');
verifyTrue(testCase, isC, 'MKAWVLK should be C-term');

% KAW: Pos 2, Prot starts with M -> N=True (Pos 2 and M)
% Len 3. End = 2+3-1 = 4. Total 7 -> C=False
[isN, isC] = srv.getWhetherProtNC('KAW');
verifyTrue(testCase, isN, 'KAW starts at 2 with M, should be N-term');
verifyFalse(testCase, isC, 'KAW ends at 4/7, not C-term');

% WVLK: Pos 4. N=False.
% Len 4. End = 4+4-1 = 7. Total 7 -> C=True
[isN, isC] = srv.getWhetherProtNC('WVLK');
verifyFalse(testCase, isN, 'WVLK starts at 4, not N-term');
verifyTrue(testCase, isC, 'WVLK ends at 7/7, should be C-term');

% XXXXX: Pos 1, Len 5, Total 5 -> N=True, C=True
[isN, isC] = srv.getWhetherProtNC('XXXXX');
verifyTrue(testCase, isN);
verifyTrue(testCase, isC);
end


function testMissingPeptide(testCase)
% TESTMISSINGPEPTIDE Validate behavior on missing peptide
% Input:
%   testCase (matlab.unittest.TestCase)
% Output:
%   (none)
srv = testCase.TestData.service;

res = srv.get_protein_name_pos('MISSING');
verifyEmpty(testCase, res);

try
    srv.getWhetherProtNC('MISSING');
    verifyFail(testCase, 'Should error for missing peptide');
catch ME
    % verify error message contains "Cannot find the peptide"
    verifyTrue(testCase, contains(ME.message, 'Cannot find the peptide'));
end
end
