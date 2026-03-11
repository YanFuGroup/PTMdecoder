classdef test_CXICPreprocessUtils < matlab.unittest.TestCase
    % Unit tests for preprocess-domain IMP group utils
    methods (Test)
        function testPreparePeakRangesBasic(testCase)
            xic_rt = (1:10)';

            range1 = struct('check_label', 1, 'rt_start', 2, 'rt_end', 4);
            current_imp_rt_range = { range1; [] };
            rt_error_tol = 0;

            [max_label, is_skip_vec, xic_peak_idx_bounds] = ...
                CXICPreprocessUtils.prepare_peak_ranges_from_imp_rt_range(...
                    xic_rt, current_imp_rt_range, rt_error_tol);

            testCase.verifyEqual(size(xic_peak_idx_bounds), [2, 1]);
            testCase.verifyEqual(size(max_label), [2, 1]);
            testCase.verifyEqual(is_skip_vec, [false; true]);
            testCase.verifyEqual(max_label(1), 1);
            testCase.verifyEqual(xic_peak_idx_bounds(1).idx_start, 2);
            testCase.verifyEqual(xic_peak_idx_bounds(1).idx_end, 4);
        end

        function testBuildPeakBoundsFromCandidatesBasic(testCase)
            xic_rt = (2:11)';
            candidate_rt_peaks = [ ...
                struct('rt_start', 2, 'rt_end', 4), ...
                struct('rt_start', 7, 'rt_end', 8)];

            [xic_peak_idx_bounds, xic_peak_rt_bounds, is_ok] = ...
                CXICPreprocessUtils.build_peak_bounds_from_candidates(...
                    xic_rt, candidate_rt_peaks);

            testCase.verifyTrue(is_ok);
            testCase.verifyEqual(size(xic_peak_idx_bounds), [2, 1]);
            testCase.verifyEqual(size(xic_peak_rt_bounds), [2, 1]);
            testCase.verifyEqual(xic_peak_idx_bounds(1).idx_start, 1);
            testCase.verifyEqual(xic_peak_idx_bounds(1).idx_end, 3);
            testCase.verifyEqual(xic_peak_idx_bounds(2).idx_start, 6);
            testCase.verifyEqual(xic_peak_idx_bounds(2).idx_end, 7);
            testCase.verifyEqual(xic_peak_rt_bounds(1).rt_start, 2);
            testCase.verifyEqual(xic_peak_rt_bounds(1).rt_end, 4);
        end

        function testBuildPeakBoundsFromCandidatesHandlesFallbackAndInvalid(testCase)
            xic_rt = (10:10:50)';
            candidate_rt_peaks = [ ...
                struct('rt_start', 8, 'rt_end', 13), ...
                struct('rt_start', 43, 'rt_end', 53), ...
                struct('rt_start', 30, 'rt_end', 20), ...
                struct('rt_start', [], 'rt_end', 40)];

            [xic_peak_idx_bounds, xic_peak_rt_bounds, is_ok] = ...
                CXICPreprocessUtils.build_peak_bounds_from_candidates(...
                    xic_rt, candidate_rt_peaks);

            testCase.verifyTrue(is_ok);
            testCase.verifyEqual(size(xic_peak_idx_bounds), [2, 1]);
            testCase.verifyEqual(xic_peak_idx_bounds(1).idx_start, 1);
            testCase.verifyEqual(xic_peak_idx_bounds(1).idx_end, 1);
            testCase.verifyEqual(xic_peak_idx_bounds(2).idx_start, 4);
            testCase.verifyEqual(xic_peak_idx_bounds(2).idx_end, 5);
            testCase.verifyEqual(xic_peak_rt_bounds(1).rt_start, 10);
            testCase.verifyEqual(xic_peak_rt_bounds(1).rt_end, 10);
            testCase.verifyEqual(xic_peak_rt_bounds(2).rt_start, 40);
            testCase.verifyEqual(xic_peak_rt_bounds(2).rt_end, 50);
        end

        function testBuildPeakBoundsFromCandidatesReturnsFalseForNoValidRanges(testCase)
            xic_rt = (1:5)';
            candidate_rt_peaks = [ ...
                struct('rt_start', 4, 'rt_end', 2), ...
                struct('rt_start', [], 'rt_end', [])];

            [xic_peak_idx_bounds, xic_peak_rt_bounds, is_ok] = ...
                CXICPreprocessUtils.build_peak_bounds_from_candidates(...
                    xic_rt, candidate_rt_peaks);

            testCase.verifyFalse(is_ok);
            testCase.verifyEmpty(xic_peak_idx_bounds);
            testCase.verifyEmpty(xic_peak_rt_bounds);
        end
    end
end
