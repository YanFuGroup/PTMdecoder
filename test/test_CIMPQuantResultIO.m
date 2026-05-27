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

            report = CIMPQuantReport();
            report = report.add_block(protein_name_pos, imp_records);
            CIMPQuantResultIO.write(report, tmp_in);

            report = CIMPQuantResultIO.read(tmp_in);
            CIMPQuantResultIO.write(report, tmp_out);

            content_in = fileread(tmp_in);
            content_out = fileread(tmp_out);
            testCase.verifyEqual(content_out, content_in);

            delete(tmp_in);
            delete(tmp_out);
        end

        function testSplitByDatasetKeepsMatchingRecordsAndHeaders(testCase)
            tmp_dir = tempname;
            mkdir(tmp_dir);
            cleanup = onCleanup(@() rmdir(tmp_dir, 's'));
            input_path = fullfile(tmp_dir, 'report_peptide_all_requant_aligned.txt');
            output_dir = fullfile(tmp_dir, 'split_by_dataset');

            lines = {
                'Protein_name,Peptide_start_position_on_protein;'
                sprintf('*\tIMP\tCharge\tDataset\tMass_center\tLow_mass_bound\tHigh_mass_bound\tPeak_area')
                sprintf('@\tRT_start\tRT_end\tProportion\tCheck_label')
                'sp|P1|PROT1,10;'
                sprintf('*\tIMP_A\t+2\tC:\\raw files\\MCF7 JIB:04.mgf\t500.0000\t499.500000\t500.500000\t100.000000')
                sprintf('@\t1.000000\t2.000000\t1.000000\t1')
                sprintf('*\tIMP_B\t+2\tMCF7_DMSO_1.mgf\t501.0000\t500.500000\t501.500000\t200.000000')
                sprintf('@\t3.000000\t4.000000\t1.000000\t0')
                'sp|P2|PROT2,20;'
                sprintf('*\tIMP_C\t+3\tMCF7_DMSO_1.mgf\t600.0000\t599.500000\t600.500000\t300.000000')
                sprintf('@\t5.000000\t6.000000\t1.000000\t1')
                };
            writeTextLines(input_path, lines);

            output_paths = CIMPQuantResultIO.splitByDataset(input_path, output_dir, 'report_peptide_all_requant_aligned');

            testCase.verifyNumElements(output_paths, 2);
            split_a = fullfile(output_dir, 'report_peptide_all_requant_aligned__MCF7_JIB_04.txt');
            split_b = fullfile(output_dir, 'report_peptide_all_requant_aligned__MCF7_DMSO_1.txt');
            testCase.verifyTrue(isfile(split_a));
            testCase.verifyTrue(isfile(split_b));

            content_a = fileread(split_a);
            testCase.verifyTrue(contains(content_a, 'Protein_name,Peptide_start_position_on_protein;'));
            testCase.verifyTrue(contains(content_a, sprintf('*\tIMP\tCharge\tDataset')));
            testCase.verifyTrue(contains(content_a, sprintf('@\tRT_start\tRT_end\tProportion\tCheck_label')));
            testCase.verifyTrue(contains(content_a, 'sp|P1|PROT1,10;'));
            testCase.verifyTrue(contains(content_a, sprintf('*\tIMP_A\t+2\tC:\\raw files\\MCF7 JIB:04.mgf')));
            testCase.verifyTrue(contains(content_a, sprintf('@\t1.000000\t2.000000\t1.000000\t1')));
            testCase.verifyFalse(contains(content_a, 'IMP_B'));
            testCase.verifyFalse(contains(content_a, 'IMP_C'));
            testCase.verifyFalse(contains(content_a, 'sp|P2|PROT2,20;'));
            testCase.verifyFalse(contains(content_a, sprintf('@\t3.000000\t4.000000\t1.000000\t0')));

            content_b = fileread(split_b);
            testCase.verifyTrue(contains(content_b, 'sp|P1|PROT1,10;'));
            testCase.verifyTrue(contains(content_b, 'sp|P2|PROT2,20;'));
            testCase.verifyTrue(contains(content_b, 'IMP_B'));
            testCase.verifyTrue(contains(content_b, 'IMP_C'));
            testCase.verifyTrue(contains(content_b, sprintf('@\t3.000000\t4.000000\t1.000000\t0')));
            testCase.verifyTrue(contains(content_b, sprintf('@\t5.000000\t6.000000\t1.000000\t1')));
            testCase.verifyFalse(contains(content_b, 'IMP_A'));

            clear cleanup;
        end
    end
end


function writeTextLines(path, lines)
% Write test report lines with newline separators.
fid = fopen(path, 'w');
cleanup = onCleanup(@() fclose(fid));
for idx_line = 1:numel(lines)
    fprintf(fid, '%s\n', lines{idx_line});
end
clear cleanup;
end
