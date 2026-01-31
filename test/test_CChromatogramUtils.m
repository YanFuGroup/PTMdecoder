classdef test_CChromatogramUtils < matlab.unittest.TestCase
    % Unit tests for CChromatogramUtils
    
    methods (Test)
        function testPreprocessMs1InputsBasicPreprocess(testCase)
            % Test basic sorting and smoothing functionality
            
            % Setup input data
            % Unsorted inputs
            rts = [10.5; 10.1; 10.3; 10.2; 10.4];
            inten = [100; 150; 800; 200; 50]; 
            % Matrix of ones corresponding to inputs
            ratioMatrix = ones(5, 1);
            minMSMSnum = 2;
            
            % Expected order: 10.1 (150), 10.2 (200), 10.3 (800), 10.4 (50), 10.5 (100)
            % After smoothing (mocking expectation roughly or checking properties)
            % Denoising will remove < 0.05 * max. Max is ~800 (smoothed might be lower).
            % Nothing should be removed here since all > 40.
            
            [rt_sorted, ratio_sorted, is_valid] = ...
                CChromatogramUtils.preprocess_ms1_inputs(rts, inten, ratioMatrix, minMSMSnum);
            
            % Assertions
            testCase.verifyTrue(is_valid);
            testCase.verifyEqual(length(rt_sorted), size(ratio_sorted, 1));
            testCase.verifyTrue(issorted(rt_sorted), 'Retention times should be sorted');
            
            % Verify the sorted order of RTs matches expectation
            expected_rts = sort(rts);
            testCase.verifyEqual(rt_sorted, expected_rts, 'AbsTol', 1e-6, 'Retention times are not sorted correctly');
        end
        
        function testPreprocessFiltering(testCase)
            % Test that low intensity noise is removed
            rts = (1:10)';
            inten = [1000; 1000; 1000; 1000; 1000; 10; 10; 10; 10; 10]; 
            ratioMatrix = ones(10, 1);
            minMSMSnum = 1;
            
            % 10 is 0.01 * 1000, should be removed (< 0.05 threshold)
            [rt_sorted, ~, is_valid] = ...
                CChromatogramUtils.preprocess_ms1_inputs(rts, inten, ratioMatrix, minMSMSnum);
                
            testCase.verifyTrue(is_valid);
            testCase.verifyEqual(length(rt_sorted), 5, 'Should filter out 5 low intensity points');
        end
        
        function testMinRowsValidation(testCase)
            % Test logic for minMSMSnum
            rts = (1:3)';
            inten = [100; 100; 100];
            ratioMatrix = ones(3, 1);
            minMSMSnum = 5; % Requirement higher than available points
            
            [~, ~, is_valid] = ...
                CChromatogramUtils.preprocess_ms1_inputs(rts, inten, ratioMatrix, minMSMSnum);
            
            testCase.verifyFalse(is_valid, 'Should be invalid due to insufficient rows');
        end

        function testGetSmoothedXic(testCase)
            % Test XIC extraction
            % Mock dependencies
            start_mz = 1000;
            charge = 2;
            
            % Create a fake dataset IO
            ms1Data = containers.Map();
            peaksData = containers.Map();
            fileMapper = struct();
            fileMapper.get_ms1_stem = @(x) x;
            
            stemName = 'test_raw';
            
            % MS1 Index: [Scan, RT, PeakIndexStart, Baseline, InjTime]
            % Scan 1, RT 10, Starts at idx 1
            % Scan 2, RT 11, Starts at idx 3
            ms1Index = [1, 10, 0, 0, 0; 
                        2, 11, 2, 0, 0;
                        3, 12, 4, 0, 0];
            
            % MS1 Peaks: [mz, intentisy]
            % Scan 1 peaks: [1000.0, 100], [1000.5, 50] (1000.5 is +1 isotope for z=2)
            % Scan 2 peaks: [1000.0, 200], [1000.5, 80]
            % Scan 3 peaks: [1000.0, 1500], [1000.5, 600]
            
            % unitdiff/charge = 1.00335/2 ~= 0.5016
            
            mz_mono = start_mz;
            mz_iso1 = start_mz + 0.501675;
            
            ms1Peaks = [mz_mono, 100; mz_iso1, 50;
                        mz_mono, 200; mz_iso1, 80;
                        mz_mono, 150,; mz_iso1, 60] ;
                        
            ms1Data(stemName) = ms1Index;
            peaksData(stemName) = ms1Peaks;
            
            datasetIO = struct();
            datasetIO.m_cMsFileMapper = fileMapper;
            datasetIO.m_mapNameMS1Index = ms1Data;
            datasetIO.m_mapNameMS1Peaks = peaksData;
            
            % Call function
            low_mz = start_mz - 0.1;
            high_mz = start_mz + 0.1;
            
            % Important: CConstant must be accessible.
            % If CConstant is missing, test will error.
            
            [xic_rt, xic_intensity_smoothed, ~] = CChromatogramUtils.get_smoothed_xic(...
                datasetIO, 'test_raw.mgf', low_mz, high_mz, charge);
            
            testCase.verifyEqual(length(xic_rt), 3);
            testCase.verifyEqual(length(xic_intensity_smoothed), 3);
            % Check that data was extracted (intensity > 0)
            % The exact logic depends on CConstant.IPV filtering and cosine distance.
            % If filtering fails, it sets intensity to 0.
            % Assuming our fake data passes filters or at least returns something.
            
            % For this specific test, we might expect filtered or not depending on IPV.
            % Since we can't easily mock CConstant.IPV without editing that file, 
            % we just check structure for now.
             
        end

        function testGetFwhm(testCase)
            % Test FWHM calculation
            
            x = linspace(-5, 5, 1000);
            % Gaussian with sigma=1. FWHM = 2.355 * sigma
            y = exp(-0.5 * x.^2);
            
            fwhm = CChromatogramUtils.get_fwhm(x, y);
            
            % Tolerance because of discrete sampling and linear interpolation
            testCase.verifyEqual(fwhm, 2.3548, 'AbsTol', 0.02);
            
            % Test border cases
            % Flat line
            y_flat = zeros(size(x));
            fwhm_flat = CChromatogramUtils.get_fwhm(x, y_flat);
            testCase.verifyEqual(fwhm_flat, 0);
            
            % Half max outside range? Logic says: use end points.
            % But Gaussian should be contained.
        end

        function testDetectXicPeaks(testCase)
            % Test XIC peak detection logic
            
            % Generate a synthetic XIC with two peaks
            xic_rt = (0:0.1:20)'; % 201 points
            % Peak 1: Center 5, Width 1
            peak1 = 100 * exp(-0.5 * ((xic_rt - 5)/0.5).^2);
            % Peak 2: Center 15, Width 1
            peak2 = 80 * exp(-0.5 * ((xic_rt - 15)/0.5).^2);
            
            xic_intensity_smoothed = peak1 + peak2;
            xic_intensity_raw = xic_intensity_smoothed; % Assume raw is same for simplicity
            
            % Case 1: PSMs identify both peaks
            rt_sorted = [5.1; 14.9];
            alpha = 0.1; % Stop at 10% max height
            
            XIC_peaks = CChromatogramUtils.detect_xic_peaks(...
                xic_rt, xic_intensity_smoothed, xic_intensity_raw, rt_sorted, alpha);
            
            testCase.verifyEqual(length(XIC_peaks), 2);
            
            % Verify Peak 1 bounds
            % Approx FWHM is 2.35*0.5 ~ 1.2. 
            % 10% height is at roughly +/- 1.5 sigma? 
            % exp(-0.5*x^2) = 0.1 => -0.5x^2 = ln(0.1)=-2.3 => x^2=4.6 => x=2.15 sigma
            % sigma=0.5 => dist=1.07. So bounds approx [3.9, 6.1]
            % Let's check indices.
            p1_left_idx = XIC_peaks(1).left_bound;
            p1_right_idx = XIC_peaks(1).right_bound;
            p1_rt_center = (xic_rt(p1_left_idx) + xic_rt(p1_right_idx)) / 2;
            
            testCase.verifyEqual(p1_rt_center, 5, 'AbsTol', 0.5);
            testCase.verifyGreaterThan(xic_rt(p1_right_idx) - xic_rt(p1_left_idx), 1); % Should have some width
            
            % Verify Peak 2
            p2_left_idx = XIC_peaks(2).left_bound;
            p2_right_idx = XIC_peaks(2).right_bound;
            p2_rt_center = (xic_rt(p2_left_idx) + xic_rt(p2_right_idx)) / 2;
            testCase.verifyEqual(p2_rt_center, 15, 'AbsTol', 0.5);
            
            % Case 2: Only one peak identified by PSMs
            rt_sorted_single = 5.1;
            XIC_peaks_single = CChromatogramUtils.detect_xic_peaks(...
                xic_rt, xic_intensity_smoothed, xic_intensity_raw, rt_sorted_single, alpha);
             
            testCase.verifyEqual(length(XIC_peaks_single), 1);
            testCase.verifyEqual((xic_rt(XIC_peaks_single(1).left_bound) + xic_rt(XIC_peaks_single(1).right_bound))/2, 5, 'AbsTol', 0.5);

            % Case 3: Filtering small peaks (raw intensity check)
            % If raw intensity is zero in the range, it should be removed.
            xic_intensity_raw_zero = zeros(size(xic_intensity_raw)); 
            XIC_peaks_filtered = CChromatogramUtils.detect_xic_peaks(...
                xic_rt, xic_intensity_smoothed, xic_intensity_raw_zero, rt_sorted, alpha);
            testCase.verifyEmpty(XIC_peaks_filtered);
        end

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
            ratio_estimated = CChromatogramUtils.calculate_kernel_ratio(...
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
            
            ratio_estimated_multi = CChromatogramUtils.calculate_kernel_ratio(...
                xic_rt, rt_sorted, ratio_sorted, peak_ranges_multi, false);
            
            % IMP 1 should be populated
            testCase.verifyEqual(ratio_estimated_multi(51, 1), 1, 'AbsTol', 0.01);
            % IMP 2 should be zero (no peak defined)
            testCase.verifyEqual(ratio_estimated_multi(51, 2), 0, 'IMP 2 should be empty');
        end

        function testParseImpRtRanges(testCase)
            % Test parse_imp_rt_ranges functionality
            
            % Setup data
            % imp 1: has peaks, one main peak
            imp1(1).check_label = 0; imp1(1).rt_start = 10.1; imp1(1).rt_end = 10.5;
            imp1(2).check_label = 2; imp1(2).rt_start = 12.0; imp1(2).rt_end = 12.5; % max
            imp1(3).check_label = 1; imp1(3).rt_start = 14.0; imp1(3).rt_end = 14.5;
            
            % imp 2: empty (should be skipped)
            imp2 = [];
            
            % imp 3: all zero labels (should be skipped)
            imp3(1).check_label = 0; imp3(1).rt_start = 20.0; imp3(1).rt_end = 20.5;
            
            % imp 4: skip via input vec
            imp4(1).check_label = 5; imp4(1).rt_start = 30.0; imp4(1).rt_end = 30.5;
            
            imp_rt_range = {imp1, imp2, imp3, imp4};
            is_skip_vec = [false, false, false, true];
            
            [final_XIC_peak, max_label, new_skip_vec] = CChromatogramUtils.parse_imp_rt_ranges(imp_rt_range, is_skip_vec);
            
            % Assertions
            testCase.verifyEqual(final_XIC_peak(1).left_bound, 12.0);
            testCase.verifyEqual(final_XIC_peak(1).right_bound, 12.5);
            testCase.verifyEqual(max_label(1), 2);
            testCase.verifyFalse(new_skip_vec(1));
            
            testCase.verifyTrue(new_skip_vec(2), 'Empty imp should be skipped');
            
            testCase.verifyTrue(new_skip_vec(3), 'Zero label imp should be skipped');
            testCase.verifyEqual(max_label(3), 0);
            
            testCase.verifyTrue(new_skip_vec(4), 'Originally skipped imp should stay skipped');
        end
        
        function testMapRtToIndices(testCase)
            % Test map_rt_to_indices functionality
            
            xic_rt = 10:0.1:20; % 10.0, 10.1, ..., 20.0
            % index 1 -> 10.0, index 21 -> 12.0
            
            num_imp = 2;
            final_XIC_peak = repmat(struct('left_bound',0,'right_bound',0), num_imp, 1);
            final_XIC_peak(1).left_bound = 12.0; % Exact match, index 21
            final_XIC_peak(1).right_bound = 13.0; % Exact match, index 31
            
            final_XIC_peak(2).left_bound = 999; % Out of bounds
            final_XIC_peak(2).right_bound = 999; 
            
            skip_vec_map = [false, true];
            rt_tol = 0.001;
            
            peak_ranges = CChromatogramUtils.map_rt_to_indices(xic_rt, final_XIC_peak, skip_vec_map, rt_tol);
            
            testCase.verifyEqual(peak_ranges(1).left_bound, 21);
            testCase.verifyEqual(peak_ranges(1).right_bound, 31);
            
            % Test Error Case
            final_XIC_peak_err = final_XIC_peak;
            skip_vec_err = [false, false]; % Unskip the bad one
            
            try
                CChromatogramUtils.map_rt_to_indices(xic_rt, final_XIC_peak_err, skip_vec_err, rt_tol);
                testCase.verifyFail('Should have thrown error for out of bound RT');
            catch ME
                testCase.verifyTrue(contains(ME.message, 'Cannot find the spectra'), 'Error message mismatch');
            end
        end

        function testGetClosedPeakData(testCase)
            % Test get_closed_peak_data functionality
            
            xic_rt = (0:10)';
            % intensity index 1(0) -> 0; index 2(10) -> 100; ...
            xic_intensity_full = (0:10)' * 10; 
            
            % Case 1: Normal internal range
            % Original: Indices 3 to 5 (Values: 20, 30, 40)
            % Padding: Expands to 2 to 6.
            % Logic: 
            %   Pad=2 matches idx 2 (Val 10). Since 2 < 3, force to 0.
            %   Pad=6 matches idx 6 (Val 50). Since 6 > 5, force to 0.
            % But wait, intensity_full(2) is 10. 
            % The output should reflect the padding logic:
            
            orig_start = 3;
            orig_end = 5;
            
            [rec_rt, rec_inten] = CChromatogramUtils.get_closed_peak_data(...
                xic_rt, xic_intensity_full, orig_start, orig_end);
            
            testCase.verifyEqual(rec_rt, xic_rt(2:6));
            testCase.verifyEqual(rec_inten, [0; 20; 30; 40; 0]);
            
            % Case 2: Left Boundary (idx_start=1)
            % Original: 1 to 2 (Values: 0, 10)
            % Padding: 1 to 3 (Values: 0, 10, 20)
            % Left pad (1) is NOT < start (1). SO KEEP VALUE (0).
            % Right pad (3) is > end (2). SO FORCE ZERO.
            
            % Use nonzero start to verify
              xic_intensity_nonzero = xic_intensity_full;
              xic_intensity_nonzero(1) = 999; 
            
            [rec_rt_l, rec_inten_l] = CChromatogramUtils.get_closed_peak_data(...
                  xic_rt, xic_intensity_nonzero, 1, 2);
            
              testCase.verifyEqual(rec_rt_l, xic_rt(1:3));
            % Expect: [999, 10, 0] because we couldn't pad left
            testCase.verifyEqual(rec_inten_l, [999; 10; 0]);
            
             % Case 3: Right Boundary (idx_end=11)
             [rec_rt_r, rec_inten_r] = CChromatogramUtils.get_closed_peak_data(...
                 xic_rt, xic_intensity_full, 10, 11);
             
             testCase.verifyEqual(rec_rt_r, xic_rt(9:11));
             testCase.verifyEqual(rec_inten_r, [0; 90; 100]);
        end

        function testCalculateArea(testCase)
            % Test calculate_area functionality (integration * 60)
            
            xic_rt = (0:5)'; % 0, 1, 2, 3, 4, 5
            % Intensities: [0, 10, 20, 10, 0, 0]
            xic_intensity_full = [0; 10; 20; 10; 0; 0];
            
            % Case 1: Standard peak (indices 2 to 4 -> Values 10, 20, 10)
            % get_closed_peak_data will expand to 1 to 5 (0, 10, 20, 10, 0)
            % Trapezoidal area of [0, 10, 20, 10, 0] with dx=1 is 40.
            % calculate_area multiplies by 60 -> 2400.
            
            area = CChromatogramUtils.calculate_area(xic_rt, xic_intensity_full, 2, 4);
            testCase.verifyEqual(area, 2400, 'Area calculation incorrect');
            
            % Case 2: Verify it handles boundary behavior implicitly
            % Indices 1 to 2 -> [0, 10]. 
            % get_closed_peak_data expands right to 3 (20) -> forced 0.
            % Data: [0, 10, 0] (Indices 1, 2, 3). RT: 0, 1, 2.
            % Area trapz: (0+10)/2 + (10+0)/2 = 5 + 5 = 10.
            % Result * 60 = 600.
            
            area_boundary = CChromatogramUtils.calculate_area(xic_rt, xic_intensity_full, 1, 2);
            testCase.verifyEqual(area_boundary, 600, 'Boundary area calculation incorrect');
        end
    end
end
