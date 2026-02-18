classdef test_CXICAlignPeakSelector < matlab.unittest.TestCase
    % Unit tests for CXICAlignPeakSelector

    methods (Test)
        function testSelectByAlignmentSizeMismatch(testCase)
            imp_max_props = ones(2, 2);
            area_imp_by_peak = ones(2, 2);
            xic_peak_rt_bounds = repmat(struct('rt_start', 0, 'rt_end', 0), 2, 2);
            rt_pred = 10; % wrong size

            testCase.verifyError(@() CXICAlignPeakSelector.select_by_alignment(...
                imp_max_props, area_imp_by_peak, xic_peak_rt_bounds, rt_pred, struct()), ...
                'select_by_alignment:InvalidRtPredSize');
        end

        function testSelectByAlignmentUsesRtWeights(testCase)
            imp_max_props = ones(2, 2);
            area_imp_by_peak = ones(2, 2);

            xic_peak_rt_bounds = repmat(struct('rt_start', 0, 'rt_end', 0), 2, 2);
            xic_peak_rt_bounds(:, 1) = struct('rt_start', 9, 'rt_end', 11);
            xic_peak_rt_bounds(:, 2) = struct('rt_start', 19, 'rt_end', 21);

            rt_pred = [10; 20];
            options = struct('rt_sigma', 0.2);

            idx_selected = CXICAlignPeakSelector.select_by_alignment(...
                imp_max_props, area_imp_by_peak, xic_peak_rt_bounds, rt_pred, options);

            testCase.verifyEqual(idx_selected, [1; 2]);
        end

        function testSelectByAlignmentIgnoresNaN(testCase)
            imp_max_props = [0.1, 0.9; 1, 1];
            area_imp_by_peak = ones(2, 2);

            xic_peak_rt_bounds = repmat(struct('rt_start', 0, 'rt_end', 0), 2, 2);
            xic_peak_rt_bounds(:, 1) = struct('rt_start', 9, 'rt_end', 11);
            xic_peak_rt_bounds(:, 2) = struct('rt_start', 19, 'rt_end', 21);

            rt_pred = [NaN; 20];
            options = struct('rt_sigma', 0.2);

            idx_selected = CXICAlignPeakSelector.select_by_alignment(...
                imp_max_props, area_imp_by_peak, xic_peak_rt_bounds, rt_pred, options);

            testCase.verifyEqual(idx_selected, [2; 2]);
        end
    end
end
