classdef test_CIMPQuantResultIO < matlab.unittest.TestCase
    % Unit tests for CIMPQuantResultIO

    methods (Test)
        function testReadWriteRoundTrip(testCase)
            tmp_in = [tempname, '.txt'];
            tmp_out = [tempname, '.txt'];

            protein_name_pos = {'P1', 10; 'P2', -1};
            rt_peaks = [
                struct('rt_start', 1.0, 'rt_end', 2.0, 'ratio', 0.5, 'check_label', 1),
                struct('rt_start', 3.0, 'rt_end', 4.0, 'ratio', 0.5, 'check_label', 0)
            ];
            imp_records = CIMPQuantRecord('IMP_A', 2, 'RAW1', 500.1234, 499.5, 500.8, 123.4, rt_peaks);

            CIMPGatherWriter.start_new_run(tmp_in);
            CIMPGatherWriter.write_imp_group_block(tmp_in, protein_name_pos, imp_records);

            report = CIMPQuantResultIO.read(tmp_in);
            CIMPQuantResultIO.write(report, tmp_out);

            content_in = fileread(tmp_in);
            content_out = fileread(tmp_out);
            testCase.verifyEqual(content_out, content_in);

            delete(tmp_in);
            delete(tmp_out);
        end
    end
end
