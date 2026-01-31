classdef test_CQuantIMPGroupAreaUtils < matlab.unittest.TestCase
    % Unit tests for area-domain IMP group utils
    methods (Test)
        function testComputeFinalArea(testCase)
            xic_rt = (1:6)';
            xic_intensity_smoothed = ones(6, 1);
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

            area_imp_final = CQuantIMPGroupAreaUtils.compute_final_area(...
                xic_rt, xic_intensity_smoothed, ratio_estimated, XIC_peaks, idx_selected);

            % Expected: imp1 uses peak1 -> 1.6*60=96, imp2 uses peak2 -> 1.2*60=72
            testCase.verifyEqual(area_imp_final, [96; 72], 'AbsTol', 1e-10);
        end

        function testGetFinalAreaFromPeakAreas(testCase)
            area_imp_by_peak = [10, 20, 30;
                                  5,  15, 25];
            idx_selected = [3; 2];

            area_imp_final = CQuantIMPGroupAreaUtils.get_final_area_from_peak_areas(...
                area_imp_by_peak, idx_selected);

            testCase.verifyEqual(area_imp_final, [30; 15], 'AbsTol', 1e-10);
        end

        function testFilterNonzeroXic(testCase)
            area_imp_final = [0; 5; 0; 2];
            rt_bound = repmat(struct('start', 0, 'end', 0), 4, 1);
            for i = 1:4
                rt_bound(i).start = i;
                rt_bound(i).end = i + 0.5;
            end
            extra_numeric = (11:14)';
            extra_empty = [];

            [imp_idx_nonzero, area_imp_final_f, rt_bound_f, extra_numeric_f, extra_empty_f] = ...
                CQuantIMPGroupAreaUtils.filter_nonzero_xic(area_imp_final, rt_bound, extra_numeric, extra_empty);

            testCase.verifyEqual(imp_idx_nonzero, [2; 4]);
            testCase.verifyEqual(area_imp_final_f, [5; 2]);
            testCase.verifyEqual([rt_bound_f.start]', [2; 4]);
            testCase.verifyEqual([rt_bound_f.end]', [2.5; 4.5]);
            testCase.verifyEqual(extra_numeric_f, [12; 14]);
            testCase.verifyEmpty(extra_empty_f);
        end

        function testFilterNonzeroXicRtOnly(testCase)
            area_imp_final = [0; 3; 4];
            rt_bound = repmat(struct('start', 0, 'end', 0), 3, 1);
            for i = 1:3
                rt_bound(i).start = i;
                rt_bound(i).end = i + 0.25;
            end

            [imp_idx_nonzero, area_imp_final_f, rt_bound_f] = ...
                CQuantIMPGroupAreaUtils.filter_nonzero_xic(area_imp_final, rt_bound);

            testCase.verifyEqual(imp_idx_nonzero, [2; 3]);
            testCase.verifyEqual(area_imp_final_f, [3; 4]);
            testCase.verifyEqual([rt_bound_f.start]', [2; 3]);
            testCase.verifyEqual([rt_bound_f.end]', [2.25; 3.25]);
        end

        function testComputeImpPeakAreaAndRatio(testCase)
            xic_rt = (1:7)';
            xic_intensity_smoothed = ones(7, 1);

            ratio_estimated = [0, 0;
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

            [area_imp_final, rt_bound, ratio_each_XIC_peak] = ...
                CQuantIMPGroupAreaUtils.compute_imp_peak_area_and_ratio(...
                    xic_rt, xic_intensity_smoothed, ratio_estimated, ...
                    peak_ranges, final_XIC_peak_for_IMP, is_skip_vec);

            testCase.verifyEqual(area_imp_final, [135; 45], 'AbsTol', 1e-10);
            testCase.verifyEqual(ratio_each_XIC_peak, [0.75; 0.25], 'AbsTol', 1e-10);
            testCase.verifyEqual([rt_bound.start]', [2; 2]);
            testCase.verifyEqual([rt_bound.end]', [4; 4]);
        end

        function testComputeImpPeakAreaAndRatioWithSkip(testCase)
            xic_rt = (1:7)';
            xic_intensity_smoothed = ones(7, 1);

            ratio_estimated = [0, 0;
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

            [area_imp_final, rt_bound, ratio_each_XIC_peak] = ...
                CQuantIMPGroupAreaUtils.compute_imp_peak_area_and_ratio(...
                    xic_rt, xic_intensity_smoothed, ratio_estimated, ...
                    peak_ranges, final_XIC_peak_for_IMP, is_skip_vec);

            testCase.verifyEqual(area_imp_final, [135; 0], 'AbsTol', 1e-10);
            testCase.verifyEqual(ratio_each_XIC_peak, [0.75; 0], 'AbsTol', 1e-10);
            testCase.verifyEqual([rt_bound.start]', [2; 0]);
            testCase.verifyEqual([rt_bound.end]', [4; 0]);
        end

        function testBuildRicFromPeaks(testCase)
            xic_rt = (1:7)';
            xic_intensity_smoothed = ones(7, 1);

            ratio_estimated = [0, 0;
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

            ric = CQuantIMPGroupAreaUtils.build_ric_from_peaks(...
                xic_rt, xic_intensity_smoothed, ratio_estimated, peak_ranges, is_skip_vec);

            testCase.verifyEqual(ric{1, 1}, xic_rt(1:5));
            testCase.verifyEqual(ric{1, 2}, [0; 1; 1; 1; 0], 'AbsTol', 1e-10);
            testCase.verifyEmpty(ric{2, 1});
            testCase.verifyEmpty(ric{2, 2});
        end
    end
end
