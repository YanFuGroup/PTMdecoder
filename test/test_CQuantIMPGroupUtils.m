classdef test_CQuantIMPGroupUtils < matlab.unittest.TestCase
    methods (Test)
        function testFilterAndNormalizePeakRatiosRemovesSmallImp(testCase)
            rt_grid = (1:5)';
            smoothed_intensity = ones(5,1);
            
            % Two IMPs: imp1 dominates, imp2 is minor in the peak range
            esti_ratio = [0, 0;
                          0.9, 0.1;
                          0.9, 0.1;
                          0.9, 0.1;
                          0, 0];
            
            XIC_peaks = struct('left_bound', 2, 'right_bound', 4);
            resFilterThres = 0.2; % imp2 should be removed (0.1/0.9 < 0.2)
            
            out_ratio = CQuantIMPGroupUtils.filter_and_normalize_peak_ratios(...
                rt_grid, smoothed_intensity, esti_ratio, XIC_peaks, resFilterThres);
            
            % Within peak range, imp2 should be zeroed and rows normalized to [1,0]
            expected_peak = [1, 0;
                             1, 0;
                             1, 0];
            testCase.verifyEqual(out_ratio(2:4, :), expected_peak, 'AbsTol', 1e-10);
            
            % Outside peak range should remain unchanged
            testCase.verifyEqual(out_ratio(1, :), [0, 0]);
            testCase.verifyEqual(out_ratio(5, :), [0, 0]);
        end
    end
end
