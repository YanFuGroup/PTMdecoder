classdef test_CIMPQuantPreprocessUtils < matlab.unittest.TestCase
    % Unit tests for preprocess-domain IMP group utils
    methods (Test)
        function testPreparePeakRangesBasic(testCase)
            xic_rt = (1:10)';

            range1 = struct('check_label', 1, 'rt_start', 2, 'rt_end', 4);
            current_imp_rt_range = { range1; [] };
            rt_error_tol = 0;

            [final_xic_peak_rt_bounds, max_label, is_skip_vec, xic_peak_idx_bounds] = ...
                CIMPQuantPreprocessUtils.prepare_peak_ranges_from_imp_rt_range(...
                    xic_rt, current_imp_rt_range, rt_error_tol);

            testCase.verifyEqual(size(final_xic_peak_rt_bounds), [2, 1]);
            testCase.verifyEqual(size(xic_peak_idx_bounds), [2, 1]);
            testCase.verifyEqual(size(max_label), [2, 1]);
            testCase.verifyEqual(is_skip_vec, [false; true]);
            testCase.verifyEqual(max_label(1), 1);
            testCase.verifyEqual(final_xic_peak_rt_bounds(1).rt_start, 2);
            testCase.verifyEqual(final_xic_peak_rt_bounds(1).rt_end, 4);
        end
    end
end
