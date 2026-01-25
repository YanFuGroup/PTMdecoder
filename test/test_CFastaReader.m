function tests = test_CFastaReader
tests = functiontests(localfunctions);
end


function testReadFasta(testCase)
% Prepare a temporary FASTA file
fastaContent = ['>sp|P12345|ProtA Description' newline ...
                'MKAWVLK' newline ...
                '>sp|P67890|ProtB Description' newline ...
                'LLQWERTY' newline ... 
                'LLQWERTY' newline ... 
                newline ...
                '>sp|P11111|ProtC' newline ...
                'AAA AA'];

filename = 'temp_test_fasta.fasta';
fid = fopen(filename, 'w');
fprintf(fid, '%s', fastaContent);
fclose(fid);

% Clean up after test
testCase.addTeardown(@() delete(filename));

% Standard regex for Uniprot: '>sp\|([^\|]+)\|'
% Or simpler if the reader uses the regex to extract the key.
% Let's assume the regex extracts the protein ID.
regex = '>sp\|([^\|]+)\|';

% Initialize Reader
reader = CFastaReader(filename, regex);

% Execute
mapProt = reader.read();

% Verify
verifyTrue(testCase, isa(mapProt, 'containers.Map'));
verifyEqual(testCase, mapProt.Count, uint64(3));

% Check content
% Note: Some readers might concat multi-line sequences or remove spaces.
% Assuming simple behavior first.
verifyTrue(testCase, isKey(mapProt, 'P12345'));
verifyEqual(testCase, mapProt('P12345'), 'MKAWVLK');

verifyTrue(testCase, isKey(mapProt, 'P67890'));
verifyEqual(testCase, mapProt('P67890'), 'LLQWERTYLLQWERTY');

verifyTrue(testCase, isKey(mapProt, 'P11111'));
verifyEqual(testCase, mapProt('P11111'), 'AAAAA');
end


function testFileNotFound(testCase)
reader = CFastaReader('non_existent.fasta', '.*');
try
    reader.read();
    verifyFail(testCase, 'Expected error for missing file was not thrown');
catch ME
    % Expected behavior if file doesn't exist
    verifyTrue(testCase, contains(ME.message, 'Can not open file'));
end
end
