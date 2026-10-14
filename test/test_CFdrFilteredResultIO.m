function tests = test_CFdrFilteredResultIO
% Regression tests for the 14-column tab-delimited FDR table contract.
tests = functiontests(localfunctions);
end

function testReadsPlainRecord(testCase)
verifyRecord(testCase, plainFields());
end

function testPreservesSpacesInIdentificationFields(testCase)
fields = plainFields();
fields{1} = 'C:\raw files\run  1\search.dat';
fields{2} = 'run  1_HCDFT.mgf';
fields{4} = 'run  1.100.100.2.dta';
fields{11} = 'protein description with spaces';
verifyRecord(testCase, fields);
end

function testPreservesModificationAndFollowingFields(testCase)
fields = plainFields();
fields{12} = 'Acetyl (Protein N-term),Phospho (STY)';
fields{13} = '0,3';
fields{14} = '112.360000';
verifyRecord(testCase, fields);
end

function testPreservesEmptyFirstField(testCase)
fields = plainFields();
fields{1} = '';
verifyRecord(testCase, fields);
end

function testPreservesConsecutiveEmptyMiddleFields(testCase)
fields = plainFields();
fields(11:13) = {'', '', ''};
verifyRecord(testCase, fields);
end

function testPreservesEmptyLastField(testCase)
fields = plainFields();
fields{14} = '';
verifyRecord(testCase, fields);
end

function testWriterReaderRoundTripPreservesTextAndEmptyFields(testCase)
fields = plainFields();
fields{1} = 'C:\raw files\run 1\search.dat';
fields{11} = '';
fields{12} = 'Phospho (STY)';
fields{14} = '';
entry = cell2struct(fields, fieldNames(), 2);
filePath = [tempname, '.txt'];
cleanup = onCleanup(@() delete(filePath)); %#ok<NASGU>
CFdrFilteredResultIO.write(entry, filePath);
result = CFdrFilteredResultIO.read(filePath);
testCase.verifyEqual(result.entries, entry);
end

function verifyRecord(testCase, fields)
filePath = [tempname, '.txt'];
cleanup = onCleanup(@() delete(filePath)); %#ok<NASGU>
fid = fopen(filePath, 'w');
assert(fid ~= -1, 'Could not create test fixture.');
closeFile = onCleanup(@() fclose(fid));
fprintf(fid, '%s\n', strjoin(fieldNames(), sprintf('\t')));
fprintf(fid, '%s\n', strjoin(fields, sprintf('\t')));
clear closeFile;
result = CFdrFilteredResultIO.read(filePath);
testCase.verifyEqual(numel(result.entries), 1);
testCase.verifyEqual(result.entries, cell2struct(fields, fieldNames(), 2));
end

function fields = plainFields()
fields = {'C:\raw\search.dat', 'run_HCDFT.mgf', '100', ...
    'run.100.100.2.dta', '2', '1000.123456', '1000.125000', ...
    '0.001544', '12', 'PEPTIDEK', 'sp|P1|PROTEIN', '-', '-', '42.000000'};
end

function names = fieldNames()
names = {'Site', 'DatasetName', 'Scan', 'Spectrum', 'Charge', ...
    'Calc_neutral_pepmass', 'precursor_neutral_mass', 'massdiff', ...
    'num_match_ions', 'peptide', 'protein', 'modification', ...
    'modificationlocation', 'Score'};
end
