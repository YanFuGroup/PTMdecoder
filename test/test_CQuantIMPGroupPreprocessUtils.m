classdef test_CQuantIMPGroupPreprocessUtils < matlab.unittest.TestCase
    % Unit tests for preprocess-domain IMP group utils
    methods (Test)
        function testPreparePeakRangesBasic(testCase)
            xic_rt = (1:10)';

            range1 = struct('check_label', 1, 'rt_start', 2, 'rt_end', 4);
            current_imp_rt_range = { range1; [] };
            rt_error_tol = 0;

            [final_XIC_peak_for_IMP, max_label, is_skip_vec, peak_ranges] = ...
                CQuantIMPGroupPreprocessUtils.prepare_peak_ranges_from_imp_rt_range(...
                    xic_rt, current_imp_rt_range, rt_error_tol);

            testCase.verifyEqual(size(final_XIC_peak_for_IMP), [2, 1]);
            testCase.verifyEqual(size(peak_ranges), [2, 1]);
            testCase.verifyEqual(size(max_label), [2, 1]);
            testCase.verifyEqual(is_skip_vec, [false; true]);
            testCase.verifyEqual(max_label(1), 1);
            testCase.verifyEqual(final_XIC_peak_for_IMP(1).left_bound, 2);
            testCase.verifyEqual(final_XIC_peak_for_IMP(1).right_bound, 4);
        end
    end
end
