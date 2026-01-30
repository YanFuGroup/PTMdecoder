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

        function testComputePeakFeaturesBasic(testCase)
            rt_grid = (1:10)';
            smoothed_intensity = ones(10, 1);
            
            % Two clearly separated peaks: [2-3] and [7-8]
            esti_ratio = [0, 0;
                          0.9, 0.1;
                          0.9, 0.1;
                          0, 0;
                          0, 0;
                          0, 0;
                          0.7, 0.3;
                          0.7, 0.3;
                          0, 0;
                          0, 0];
            
            XIC_peaks = repmat(struct('left_bound', 0, 'right_bound', 0), 1, 2);
            XIC_peaks(1).left_bound = 2;
            XIC_peaks(1).right_bound = 3;
            XIC_peaks(2).left_bound = 7;
            XIC_peaks(2).right_bound = 8;
            
            [imp_max_props, peak_fwhms, ratio_each_XIC_peak, rt_bound] = ...
                CQuantIMPGroupUtils.compute_peak_features(rt_grid, smoothed_intensity, esti_ratio, XIC_peaks);
            
            testCase.verifyEqual(size(imp_max_props), [2, 2]);
            testCase.verifyEqual(size(peak_fwhms), [2, 2]);
            testCase.verifyEqual(size(ratio_each_XIC_peak), [2, 2]);
            testCase.verifyEqual(size(rt_bound), [2, 2]);
            
            testCase.verifyEqual(imp_max_props(:, 1), [0.9; 0.1], 'AbsTol', 1e-10);
            testCase.verifyEqual(imp_max_props(:, 2), [0.7; 0.3], 'AbsTol', 1e-10);
            testCase.verifyEqual(peak_fwhms(:, 1), [1; 1], 'AbsTol', 1e-10);
            testCase.verifyEqual(peak_fwhms(:, 2), [1; 1], 'AbsTol', 1e-10);
            
            % Expected areas: peak1 imp1 -> 1.8*60=108, imp2 -> 0.2*60=12
            %                 peak2 imp1 -> 1.4*60=84,  imp2 -> 0.6*60=36
            testCase.verifyEqual(ratio_each_XIC_peak(:, 1), [108; 12], 'AbsTol', 1e-10);
            testCase.verifyEqual(ratio_each_XIC_peak(:, 2), [84; 36], 'AbsTol', 1e-10);
            
            testCase.verifyEqual(rt_bound(1, 1).start, 2);
            testCase.verifyEqual(rt_bound(1, 1).end, 3);
            testCase.verifyEqual(rt_bound(1, 2).start, 7);
            testCase.verifyEqual(rt_bound(1, 2).end, 8);
            testCase.verifyEqual(rt_bound(2, 1).start, 2);
            testCase.verifyEqual(rt_bound(2, 1).end, 3);
            testCase.verifyEqual(rt_bound(2, 2).start, 7);
            testCase.verifyEqual(rt_bound(2, 2).end, 8);
        end

        function testSelectBestPeakPerImp(testCase)
            imp_max_props = [0.9, 0.2;
                             0.1, 0.8];
            ratio_each_XIC_peak = [100, 10;
                                   5, 60];
            % Scores:
            % imp1 -> [90, 2] -> pick 1
            % imp2 -> [0.5, 48] -> pick 2
            idx_selected = CQuantIMPGroupUtils.select_best_peak_per_imp(...
                imp_max_props, ratio_each_XIC_peak);
            
            testCase.verifyEqual(idx_selected, [1; 2]);
        end
    end
end
