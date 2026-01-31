classdef test_CQuantIMPGroupUtils < matlab.unittest.TestCase
    % Unit tests for CQuantIMPGroupUtils
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
            
            [imp_max_props, peak_fwhms, area_each_XIC_peak, rt_bound] = ...
                CQuantIMPGroupUtils.compute_peak_features(rt_grid, smoothed_intensity, esti_ratio, XIC_peaks);
            
            testCase.verifyEqual(size(imp_max_props), [2, 2]);
            testCase.verifyEqual(size(peak_fwhms), [2, 2]);
            testCase.verifyEqual(size(area_each_XIC_peak), [2, 2]);
            testCase.verifyEqual(size(rt_bound), [2, 2]);
            
            testCase.verifyEqual(imp_max_props(:, 1), [0.9; 0.1], 'AbsTol', 1e-10);
            testCase.verifyEqual(imp_max_props(:, 2), [0.7; 0.3], 'AbsTol', 1e-10);
            testCase.verifyEqual(peak_fwhms(:, 1), [1; 1], 'AbsTol', 1e-10);
            testCase.verifyEqual(peak_fwhms(:, 2), [1; 1], 'AbsTol', 1e-10);
            
            % Expected areas: peak1 imp1 -> 1.8*60=108, imp2 -> 0.2*60=12
            %                 peak2 imp1 -> 1.4*60=84,  imp2 -> 0.6*60=36
            testCase.verifyEqual(area_each_XIC_peak(:, 1), [108; 12], 'AbsTol', 1e-10);
            testCase.verifyEqual(area_each_XIC_peak(:, 2), [84; 36], 'AbsTol', 1e-10);
            
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
            area_each_XIC_peak = [100, 10;
                                  5, 60];
            % Scores:
            % imp1 -> [90, 2] -> pick 1
            % imp2 -> [0.5, 48] -> pick 2
            idx_selected = CQuantIMPGroupUtils.select_best_peak_per_imp(...
                imp_max_props, area_each_XIC_peak);
            
            testCase.verifyEqual(idx_selected, [1; 2]);
        end

        function testRefineRatiosBySelection(testCase)
            esti_ratio = [0, 0;
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
            
            out_ratio = CQuantIMPGroupUtils.refine_ratios_by_selection(...
                esti_ratio, XIC_peaks, idx_selected);
            
            % IMP1 keeps peak1 only (rows 2-3)
            testCase.verifyEqual(out_ratio(:, 1), [0; 0.8; 0.8; 0; 0; 0], 'AbsTol', 1e-10);
            % IMP2 keeps peak2 only (rows 4-5)
            testCase.verifyEqual(out_ratio(:, 2), [0; 0; 0; 0.6; 0.6; 0], 'AbsTol', 1e-10);
        end

        function testComputeFinalArea(testCase)
            rt_grid = (1:6)';
            smoothed_intensity = ones(6, 1);
            esti_ratio = [0, 0;
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
            
            auxic = CQuantIMPGroupUtils.compute_final_area(...
                rt_grid, smoothed_intensity, esti_ratio, XIC_peaks, idx_selected);
            
            % Expected: imp1 uses peak1 -> 1.6*60=96, imp2 uses peak2 -> 1.2*60=72
            testCase.verifyEqual(auxic, [96; 72], 'AbsTol', 1e-10);
        end

        function testGetFinalAreaFromPeakAreas(testCase)
            area_each_XIC_peak = [10, 20, 30;
                                  5,  15, 25];
            idx_selected = [3; 2];

            auxic = CQuantIMPGroupUtils.get_final_area_from_peak_areas(...
                area_each_XIC_peak, idx_selected);

            testCase.verifyEqual(auxic, [30; 15], 'AbsTol', 1e-10);
        end

        function testFilterNonzeroXic(testCase)
            auxic = [0; 5; 0; 2];
            rt_bound = repmat(struct('start', 0, 'end', 0), 4, 1);
            for i = 1:4
                rt_bound(i).start = i;
                rt_bound(i).end = i + 0.5;
            end
            extra_numeric = (11:14)';
            extra_empty = [];

            [idxNonZero, auxic_f, rt_bound_f, extra_numeric_f, extra_empty_f] = ...
                CQuantIMPGroupUtils.filter_nonzero_xic(auxic, rt_bound, extra_numeric, extra_empty);

            testCase.verifyEqual(idxNonZero, [2; 4]);
            testCase.verifyEqual(auxic_f, [5; 2]);
            testCase.verifyEqual([rt_bound_f.start]', [2; 4]);
            testCase.verifyEqual([rt_bound_f.end]', [2.5; 4.5]);
            testCase.verifyEqual(extra_numeric_f, [12; 14]);
            testCase.verifyEmpty(extra_empty_f);
        end

        function testFilterNonzeroXicRtOnly(testCase)
            auxic = [0; 3; 4];
            rt_bound = repmat(struct('start', 0, 'end', 0), 3, 1);
            for i = 1:3
                rt_bound(i).start = i;
                rt_bound(i).end = i + 0.25;
            end

            [idxNonZero, auxic_f, rt_bound_f] = ...
                CQuantIMPGroupUtils.filter_nonzero_xic(auxic, rt_bound);

            testCase.verifyEqual(idxNonZero, [2; 3]);
            testCase.verifyEqual(auxic_f, [3; 4]);
            testCase.verifyEqual([rt_bound_f.start]', [2; 3]);
            testCase.verifyEqual([rt_bound_f.end]', [2.25; 3.25]);
        end

        function testComputeImpPeakAreaAndRatio(testCase)
            rt_grid = (1:7)';
            smoothed_intensity = ones(7, 1);

            esti_ratio = [0, 0;
                          0.75, 0.25;
                          0.75, 0.25;
                          0.75, 0.25;
                          0, 0;
                          0.2, 0.8;
                          0, 0];

            peak_ranges = repmat(struct('left_bound', 0, 'right_bound', 0), 1, 2);
            peak_ranges(1).left_bound = 2;
            peak_ranges(1).right_bound = 4;
            peak_ranges(2).left_bound = 2;
            peak_ranges(2).right_bound = 4;

            final_XIC_peak_for_IMP = repmat(struct('left_bound', 0, 'right_bound', 0), 1, 2);
            final_XIC_peak_for_IMP(1).left_bound = 2;
            final_XIC_peak_for_IMP(1).right_bound = 4;
            final_XIC_peak_for_IMP(2).left_bound = 2;
            final_XIC_peak_for_IMP(2).right_bound = 4;

            is_skip_vec = [false; false];

            [auxic, rt_bound, ratio_each_XIC_peak] = ...
                CQuantIMPGroupUtils.compute_imp_peak_area_and_ratio(...
                    rt_grid, smoothed_intensity, esti_ratio, ...
                    peak_ranges, final_XIC_peak_for_IMP, is_skip_vec);

            testCase.verifyEqual(auxic, [135; 45], 'AbsTol', 1e-10);
            testCase.verifyEqual(ratio_each_XIC_peak, [0.75; 0.25], 'AbsTol', 1e-10);
            testCase.verifyEqual([rt_bound.start]', [2; 2]);
            testCase.verifyEqual([rt_bound.end]', [4; 4]);
        end

        function testComputeImpPeakAreaAndRatioWithSkip(testCase)
            rt_grid = (1:7)';
            smoothed_intensity = ones(7, 1);

            esti_ratio = [0, 0;
                          0.75, 0.25;
                          0.75, 0.25;
                          0.75, 0.25;
                          0, 0;
                          0.2, 0.8;
                          0, 0];

            peak_ranges = repmat(struct('left_bound', 0, 'right_bound', 0), 1, 2);
            peak_ranges(1).left_bound = 2;
            peak_ranges(1).right_bound = 4;
            peak_ranges(2).left_bound = 2;
            peak_ranges(2).right_bound = 4;

            final_XIC_peak_for_IMP = repmat(struct('left_bound', 0, 'right_bound', 0), 1, 2);
            final_XIC_peak_for_IMP(1).left_bound = 2;
            final_XIC_peak_for_IMP(1).right_bound = 4;
            final_XIC_peak_for_IMP(2).left_bound = 2;
            final_XIC_peak_for_IMP(2).right_bound = 4;

            is_skip_vec = [false; true];

            [auxic, rt_bound, ratio_each_XIC_peak] = ...
                CQuantIMPGroupUtils.compute_imp_peak_area_and_ratio(...
                    rt_grid, smoothed_intensity, esti_ratio, ...
                    peak_ranges, final_XIC_peak_for_IMP, is_skip_vec);

            testCase.verifyEqual(auxic, [135; 0], 'AbsTol', 1e-10);
            testCase.verifyEqual(ratio_each_XIC_peak, [0.75; 0], 'AbsTol', 1e-10);
            testCase.verifyEqual([rt_bound.start]', [2; 0]);
            testCase.verifyEqual([rt_bound.end]', [4; 0]);
        end

        function testBuildRicFromPeaks(testCase)
            rt_grid = (1:7)';
            smoothed_intensity = ones(7, 1);

            esti_ratio = [0, 0;
                          1, 0;
                          1, 0;
                          1, 0;
                          0, 0;
                          0.2, 0.8;
                          0, 0];

            peak_ranges = repmat(struct('left_bound', 0, 'right_bound', 0), 1, 2);
            peak_ranges(1).left_bound = 2;
            peak_ranges(1).right_bound = 4;
            peak_ranges(2).left_bound = 2;
            peak_ranges(2).right_bound = 4;

            is_skip_vec = [false; true];

            ric = CQuantIMPGroupUtils.build_ric_from_peaks(...
                rt_grid, smoothed_intensity, esti_ratio, peak_ranges, is_skip_vec);

            testCase.verifyEqual(ric{1, 1}, rt_grid(1:5));
            testCase.verifyEqual(ric{1, 2}, [0; 1; 1; 1; 0], 'AbsTol', 1e-10);
            testCase.verifyEmpty(ric{2, 1});
            testCase.verifyEmpty(ric{2, 2});
        end
    end
end
