classdef test_CIMPGatherWriter < matlab.unittest.TestCase
    % Unit tests for CIMPGatherWriter
    methods (Test)
        function testWriteAllSingleRecord(testCase)
            tmp_file = [tempname, '.txt'];
            protein_name_pos = {'P1', 10};
            rt_peaks = struct('rt_start', 1.0, 'rt_end', 2.0, 'ratio', 0.5, 'check_label', 1);
            imp_records = CIMPQuantResult('IMP_A', 2, 'RAW1', 500.1234, 499.5, 500.8, 123.4, rt_peaks);

            CIMPGatherWriter.start_new_run(tmp_file);
            CIMPGatherWriter.write_imp_group_block(tmp_file, protein_name_pos, imp_records);
            content = fileread(tmp_file);

            testCase.verifyTrue(contains(content, 'Protein_name,Peptide_start_position_on_protein;'));
            testCase.verifyTrue(contains(content, 'P1,10;'));
            testCase.verifyTrue(~isempty(regexp(content, '\*\tIMP_A\t\+2\tRAW1', 'once')));
            testCase.verifyTrue(~isempty(regexp(content, '@\t1\.000000\t2\.000000\t0\.500000\t1', 'once')));

            delete(tmp_file);
        end

        function testWriteAllEmptyRecordsNoFile(testCase)
            tmp_file = [tempname, '.txt'];
            protein_name_pos = {'P1', 10};
            imp_records = CIMPQuantResult.empty(0,1);

            CIMPGatherWriter.start_new_run(tmp_file);
            CIMPGatherWriter.write_imp_group_block(tmp_file, protein_name_pos, imp_records);
            testCase.verifyEqual(exist(tmp_file, 'file'), 2);
            file_info = dir(tmp_file);
            testCase.verifyGreaterThan(file_info.bytes, 0);
            delete(tmp_file);
        end
    end
end
