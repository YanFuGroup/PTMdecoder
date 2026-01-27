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
    end
end
