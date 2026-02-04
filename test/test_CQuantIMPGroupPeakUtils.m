classdef test_CQuantIMPGroupPeakUtils < matlab.unittest.TestCase
    % Unit tests for peak-domain IMP group utils
    methods (Test)
        function testCalculateKernelRatio(testCase)
            % Test Kernel Ratio Calculation

            % Grid: 0 to 10
            xic_rt = (0:0.1:10)';

            % PSMs at RT=5. Ratio=0.5 for IMP1, 0.8 for IMP2
            rt_sorted = [4.9; 5.0; 5.1];
            ratio_sorted = [0.2, 0.8;
                            0.2, 0.8;
                            0.2, 0.8];

            % Define one peak range covering these PSMs (e.g. 4.0 to 6.0)
            % Indices: 4.0 is index 41, 6.0 is index 61
            peak_range = struct('left_bound', 41, 'right_bound', 61);

            % SCENARIO 1: is_broadcast = true (Same peak logic for all IMPs)
            % Should produce ratios approx 0.5 and 0.8 in the peak region
            ratio_estimated = CQuantIMPGroupPeakUtils.calculate_kernel_ratio(...
                xic_rt, rt_sorted, ratio_sorted, peak_range, true);

            % Check dimensions
            testCase.verifyEqual(size(ratio_estimated), [101, 2]);

            % Check value at center (index 51, rt=5.0)
            center_ratio = ratio_estimated(51, :);
            % Since all PSMs have same ratio, the weighted average should be exactly that ratio.
            % (Weights are normalized)
            testCase.verifyEqual(center_ratio(1), 0.2, 'AbsTol', 0.01);
            testCase.verifyEqual(center_ratio(2), 0.8, 'AbsTol', 0.01);

            % Check that outside the peak, it is zero
            testCase.verifyEqual(ratio_estimated(10, :), [0 0]);

            % SCENARIO 2: is_broadcast = false (Different peaks per IMP)
            % IMP 1 uses the peak at 5.0. IMP 2 has NO peak (or different peak).

            peak_ranges_multi = repmat(struct('left_bound',0,'right_bound',0), 1, 2);
            peak_ranges_multi(1) = peak_range; % IMP 1 has peak
            % IMP 2 is empty/default

            ratio_estimated_multi = CQuantIMPGroupPeakUtils.calculate_kernel_ratio(...
                xic_rt, rt_sorted, ratio_sorted, peak_ranges_multi, false);

            % IMP 1 should be populated
            testCase.verifyEqual(ratio_estimated_multi(51, 1), 1, 'AbsTol', 0.01);
            % IMP 2 should be zero (no peak defined)
            testCase.verifyEqual(ratio_estimated_multi(51, 2), 0, 'IMP 2 should be empty');
        end
        function testFilterAndNormalizePeakRatiosRemovesSmallImp(testCase)
            xic_rt = (1:5)';
            xic_intensity_smoothed = ones(5,1);

            % Two IMPs: imp1 dominates, imp2 is minor in the peak range
            ratio_estimated = [0, 0;
                          0.9, 0.1;
                          0.9, 0.1;
                          0.9, 0.1;
                          0, 0];

            XIC_peaks = struct('left_bound', 2, 'right_bound', 4);
            resFilterThres = 0.2; % imp2 should be removed (0.1/0.9 < 0.2)

            ratio_estimated_out = CQuantIMPGroupPeakUtils.filter_and_normalize_peak_ratios(...
                xic_rt, xic_intensity_smoothed, ratio_estimated, XIC_peaks, resFilterThres);

            % Within peak range, imp2 should be zeroed and rows normalized to [1,0]
            expected_peak = [1, 0;
                             1, 0;
                             1, 0];
            testCase.verifyEqual(ratio_estimated_out(2:4, :), expected_peak, 'AbsTol', 1e-10);

            % Outside peak range should remain unchanged
            testCase.verifyEqual(ratio_estimated_out(1, :), [0, 0]);
            testCase.verifyEqual(ratio_estimated_out(5, :), [0, 0]);
        end

        function testComputePeakFeaturesBasic(testCase)
            xic_rt = (1:10)';
            xic_intensity_smoothed = ones(10, 1);

            % Two clearly separated peaks: [2-3] and [7-8]
            ratio_estimated = [0, 0;
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

            [imp_max_props, peak_fwhms, area_imp_by_peak, rt_bound] = ...
                CQuantIMPGroupPeakUtils.compute_peak_features(xic_rt, xic_intensity_smoothed, ratio_estimated, XIC_peaks);

            testCase.verifyEqual(size(imp_max_props), [2, 2]);
            testCase.verifyEqual(size(peak_fwhms), [2, 2]);
            testCase.verifyEqual(size(area_imp_by_peak), [2, 2]);
            testCase.verifyEqual(size(rt_bound), [2, 2]);

            testCase.verifyEqual(imp_max_props(:, 1), [0.9; 0.1], 'AbsTol', 1e-10);
            testCase.verifyEqual(imp_max_props(:, 2), [0.7; 0.3], 'AbsTol', 1e-10);
            testCase.verifyEqual(peak_fwhms(:, 1), [1; 1], 'AbsTol', 1e-10);
            testCase.verifyEqual(peak_fwhms(:, 2), [1; 1], 'AbsTol', 1e-10);

            % Expected areas: peak1 imp1 -> 1.8*60=108, imp2 -> 0.2*60=12
            %                 peak2 imp1 -> 1.4*60=84,  imp2 -> 0.6*60=36
            testCase.verifyEqual(area_imp_by_peak(:, 1), [108; 12], 'AbsTol', 1e-10);
            testCase.verifyEqual(area_imp_by_peak(:, 2), [84; 36], 'AbsTol', 1e-10);

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
            area_imp_by_peak = [100, 10;
                                  5, 60];
            % Scores:
            % imp1 -> [90, 2] -> pick 1
            % imp2 -> [0.5, 48] -> pick 2
            idx_selected = CQuantIMPGroupPeakUtils.select_best_peak_per_imp(...
                imp_max_props, area_imp_by_peak);

            testCase.verifyEqual(idx_selected, [1; 2]);
        end

        function testRefineRatiosBySelection(testCase)
            ratio_estimated = [0, 0;
                          0.8, 0.2;
                          0.8, 0.2;
                          0.4, 0.6;
                          0.4, 0.6;
                          0, 0];
            XIC_peaks = repmat(struct('left_bound', 0, 'right_bound', 0), 1, 2);
            XIC_peaks(1).left_bound = 2;
            XIC_peaks(1).right_bound = 3;
            XIC_peaks(2).left_bound = 4;
            XIC_peaks(2).right_bound = 5;
            idx_selected = [1; 2];

            ratio_estimated_out = CQuantIMPGroupPeakUtils.refine_ratios_by_selection(...
                ratio_estimated, XIC_peaks, idx_selected);

            % IMP1 keeps peak1 only (rows 2-3)
            testCase.verifyEqual(ratio_estimated_out(:, 1), [0; 0.8; 0.8; 0; 0; 0], 'AbsTol', 1e-10);
            % IMP2 keeps peak2 only (rows 4-5)
            testCase.verifyEqual(ratio_estimated_out(:, 2), [0; 0; 0; 0.6; 0.6; 0], 'AbsTol', 1e-10);
        end

        function testHasMinRows(testCase)
            % Test hasMinRows function

            % Test with default min_rows = 1
            ratio_matrix_empty = zeros(0, 2);
            testCase.verifyFalse(CQuantIMPGroupPeakUtils.hasMinRows(ratio_matrix_empty));

            ratio_matrix_one_row = ones(1, 2);
            testCase.verifyTrue(CQuantIMPGroupPeakUtils.hasMinRows(ratio_matrix_one_row));

            ratio_matrix_multi_row = ones(5, 2);
            testCase.verifyTrue(CQuantIMPGroupPeakUtils.hasMinRows(ratio_matrix_multi_row));

            % Test with specified min_rows
            testCase.verifyTrue(CQuantIMPGroupPeakUtils.hasMinRows(ratio_matrix_multi_row, 3));
            testCase.verifyFalse(CQuantIMPGroupPeakUtils.hasMinRows(ratio_matrix_one_row, 3));
            testCase.verifyTrue(CQuantIMPGroupPeakUtils.hasMinRows(ratio_matrix_multi_row, 5));
            testCase.verifyFalse(CQuantIMPGroupPeakUtils.hasMinRows(ratio_matrix_multi_row, 6));
        end
    end
end
