classdef test_CXICAlignScore < matlab.unittest.TestCase
    % Unit tests for CXICAlignScore

    methods (Test)
        function testWeightsSymmetryDefaultSigma(testCase)
            xic_peak_rt_bounds = repmat(struct('rt_start', 0, 'rt_end', 0), 1, 2);
            xic_peak_rt_bounds(1).rt_start = 4;
            xic_peak_rt_bounds(1).rt_end = 4;
            xic_peak_rt_bounds(2).rt_start = 6;
            xic_peak_rt_bounds(2).rt_end = 6;

            weights = CXICAlignScore.compute_rt_weights(5, xic_peak_rt_bounds, []);

            testCase.verifyEqual(weights(1), weights(2), 'AbsTol', 1e-12);
        end

        function testWeightsMonotonicity(testCase)
            xic_peak_rt_bounds = repmat(struct('rt_start', 0, 'rt_end', 0), 1, 2);
            xic_peak_rt_bounds(1).rt_start = 6;
            xic_peak_rt_bounds(1).rt_end = 6;
            xic_peak_rt_bounds(2).rt_start = 8;
            xic_peak_rt_bounds(2).rt_end = 8;

            weights = CXICAlignScore.compute_rt_weights(5, xic_peak_rt_bounds, 0.5);

            testCase.verifyGreaterThan(weights(1), weights(2));
        end
    end
end
