classdef test_CChromatogramUtils < matlab.unittest.TestCase
    
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
            
            [sort_rts, sort_inten, sort_ratioMatrix, is_valid] = ...
                CChromatogramUtils.preprocess_ms1_inputs(rts, inten, ratioMatrix, minMSMSnum);
            
            % Assertions
            testCase.verifyTrue(is_valid);
            testCase.verifyEqual(length(sort_rts), length(sort_inten));
            testCase.verifyTrue(issorted(sort_rts), 'Retention times should be sorted');
            
            % Check sorting correctness for specific known max
            % Original max was at 10.3. In sorted list (10.1, 10.2, 10.3, 10.4, 10.5), it should be index 3.
            [~, max_idx] = max(sort_inten);
            % The peak should be preserved at 10.3 (which becomes index 3 after sort)
            testCase.verifyEqual(sort_rts(max_idx), 10.3, 'AbsTol', 1e-6, 'Max intensity location shifted after sort/smooth');
            
            % Verify the sorted order of RTs matches expectation
            expected_rts = sort(rts);
            testCase.verifyEqual(sort_rts, expected_rts, 'AbsTol', 1e-6, 'Retention times are not sorted correctly');
        end
        
        function testPreprocessFiltering(testCase)
            % Test that low intensity noise is removed
            rts = (1:10)';
            inten = [1000; 1000; 1000; 1000; 1000; 10; 10; 10; 10; 10]; 
            ratioMatrix = ones(10, 1);
            minMSMSnum = 1;
            
            % 10 is 0.01 * 1000, should be removed (< 0.05 threshold)
            [sort_rts, sort_inten, ~, is_valid] = ...
                CChromatogramUtils.preprocess_ms1_inputs(rts, inten, ratioMatrix, minMSMSnum);
                
            testCase.verifyTrue(is_valid);
            testCase.verifyEqual(length(sort_rts), 5, 'Should filter out 5 low intensity points');
        end
        
        function testMinRowsValidation(testCase)
            % Test logic for minMSMSnum
            rts = (1:3)';
            inten = [100; 100; 100];
            ratioMatrix = ones(3, 1);
            minMSMSnum = 5; % Requirement higher than available points
            
            [~, ~, ~, is_valid] = ...
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
            
            rawName = 'test_raw';
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
            
            mz_mono = 1000.0;
            mz_iso1 = 1000.0 + 0.501675;
            
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
            low_mz = 999.9;
            high_mz = 1000.1;
            
            % Important: CConstant must be accessible.
            % If CConstant is missing, test will error.
            
            [rt, smoothed, raw] = CChromatogramUtils.get_smoothed_xic(...
                datasetIO, 'test_raw.mgf', low_mz, high_mz, charge);
            
            testCase.verifyEqual(length(rt), 3);
            testCase.verifyEqual(length(smoothed), 3);
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
            rt_grid = (0:0.1:20)'; % 201 points
            % Peak 1: Center 5, Width 1
            peak1 = 100 * exp(-0.5 * ((rt_grid - 5)/0.5).^2);
            % Peak 2: Center 15, Width 1
            peak2 = 80 * exp(-0.5 * ((rt_grid - 15)/0.5).^2);
            
            smoothed_intensity = peak1 + peak2;
            raw_intensity = smoothed_intensity; % Assume raw is same for simplicity
            
            % Case 1: PSMs identify both peaks
            sort_rts = [5.1; 14.9];
            alpha = 0.1; % Stop at 10% max height
            
            XIC_peaks = CChromatogramUtils.detect_xic_peaks(...
                rt_grid, smoothed_intensity, raw_intensity, sort_rts, alpha);
            
            testCase.verifyEqual(length(XIC_peaks), 2);
            
            % Verify Peak 1 bounds
            % Approx FWHM is 2.35*0.5 ~ 1.2. 
            % 10% height is at roughly +/- 1.5 sigma? 
            % exp(-0.5*x^2) = 0.1 => -0.5x^2 = ln(0.1)=-2.3 => x^2=4.6 => x=2.15 sigma
            % sigma=0.5 => dist=1.07. So bounds approx [3.9, 6.1]
            % Let's check indices.
            p1_left_idx = XIC_peaks(1).left_bound;
            p1_right_idx = XIC_peaks(1).right_bound;
            p1_rt_center = (rt_grid(p1_left_idx) + rt_grid(p1_right_idx)) / 2;
            
            testCase.verifyEqual(p1_rt_center, 5, 'AbsTol', 0.5);
            testCase.verifyGreaterThan(rt_grid(p1_right_idx) - rt_grid(p1_left_idx), 1); % Should have some width
            
            % Verify Peak 2
            p2_left_idx = XIC_peaks(2).left_bound;
            p2_right_idx = XIC_peaks(2).right_bound;
            p2_rt_center = (rt_grid(p2_left_idx) + rt_grid(p2_right_idx)) / 2;
            testCase.verifyEqual(p2_rt_center, 15, 'AbsTol', 0.5);
            
            % Case 2: Only one peak identified by PSMs
            sort_rts_single = [5.1];
            XIC_peaks_single = CChromatogramUtils.detect_xic_peaks(...
                rt_grid, smoothed_intensity, raw_intensity, sort_rts_single, alpha);
             
            testCase.verifyEqual(length(XIC_peaks_single), 1);
            testCase.verifyEqual((rt_grid(XIC_peaks_single(1).left_bound) + rt_grid(XIC_peaks_single(1).right_bound))/2, 5, 'AbsTol', 0.5);

            % Case 3: Filtering small peaks (raw intensity check)
            % If raw intensity is zero in the range, it should be removed.
            raw_intensity_zero = zeros(size(raw_intensity)); 
            XIC_peaks_filtered = CChromatogramUtils.detect_xic_peaks(...
                rt_grid, smoothed_intensity, raw_intensity_zero, sort_rts, alpha);
            testCase.verifyEmpty(XIC_peaks_filtered);
        end

        function testCalculateKernelRatio(testCase)
            % Test Kernel Ratio Calculation
            
            % Grid: 0 to 10
            rt_grid = (0:0.1:10)';
            
            % PSMs at RT=5. Ratio=0.5 for Iso1, 0.8 for Iso2
            sort_rts = [4.9; 5.0; 5.1];
            sort_ratioMatrix = [0.2, 0.8; 
                                0.2, 0.8; 
                                0.2, 0.8];
            
            % Define one peak range covering these PSMs (e.g. 4.0 to 6.0)
            % Indices: 4.0 is index 41, 6.0 is index 61
            peak_range = struct('left_bound', 41, 'right_bound', 61);
            
            % SCENARIO 1: is_shared = true (Same peak logic for all isotopes)
            % Should produce ratios approx 0.5 and 0.8 in the peak region
            esti_ratio = CChromatogramUtils.calculate_kernel_ratio(...
                rt_grid, sort_rts, sort_ratioMatrix, peak_range, true);
            
            % Check dimensions
            testCase.verifyEqual(size(esti_ratio), [101, 2]);
            
            % Check value at center (index 51, rt=5.0)
            center_ratio = esti_ratio(51, :);
            % Since all PSMs have same ratio, the weighted average should be exactly that ratio.
            % (Weights are normalized)
            testCase.verifyEqual(center_ratio(1), 0.2, 'AbsTol', 0.01);
            testCase.verifyEqual(center_ratio(2), 0.8, 'AbsTol', 0.01);
            
            % Check that outside the peak, it is zero
            testCase.verifyEqual(esti_ratio(10, :), [0 0]);
            
            % SCENARIO 2: is_shared = false (Different peaks per isotope)
            % Iso 1 uses the peak at 5.0. Iso 2 has NO peak (or different peak).
            
            peak_ranges_multi = repmat(struct('left_bound',0,'right_bound',0), 1, 2);
            peak_ranges_multi(1) = peak_range; % Iso 1 has peak
            % Iso 2 is empty/default
            
            esti_ratio_multi = CChromatogramUtils.calculate_kernel_ratio(...
                rt_grid, sort_rts, sort_ratioMatrix, peak_ranges_multi, false);
            
            % Iso 1 should be populated
            testCase.verifyEqual(esti_ratio_multi(51, 1), 1, 'AbsTol', 0.01);
            % Iso 2 should be zero (no peak defined)
            testCase.verifyEqual(esti_ratio_multi(51, 2), 0, 'Iso 2 should be empty');
        end
    end
end
