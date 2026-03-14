classdef test_estimateDatasetNoiseModel < matlab.unittest.TestCase
% Unit test for CMS2QuantSolver.estimateDatasetNoiseModel

methods (Test)

    function testNormalDataset(testCase)
        % Test basic pooled MLE and MAD calculation behavior
        
        % Setup mock inputs
        inputs(1).filteredOutExpPeakCount = 100;
        inputs(1).filteredOutExpPeakSqSum = 400; % variance = 4 => sigma = 2
        inputs(1).matchedExpPeaks = [1, 100, 100; 2, 80, 80]; % [ionIndex, normInt, rawInt]
        inputs(1).fittedMatchedPeakIntensities = [90; 80]; % First has residual 10, second 0
        
        inputs(2).filteredOutExpPeakCount = 50;
        inputs(2).filteredOutExpPeakSqSum = 200; % Total count=150, SqSum=600 => sigma=2
        inputs(2).matchedExpPeaks = [1, 200, 200; 2, 5, 5]; 
        inputs(2).fittedMatchedPeakIntensities = [220; 5]; % First has residual -20, second 0
        
        % Execute
        model = CMS2QuantSolver.estimateDatasetNoiseModel(inputs);
        
        % Verify sigma_base
        testCase.verifyEqual(model.sigma_base, 2.0, 'AbsTol', 1e-4);
        testCase.verifyEqual(model.tau_floor, 10.0, 'AbsTol', 1e-4);
        
        % Theoretical high signals: 90, 80, 220 (all > 10)
        % Raw y: 100, 80, 200
        % Residuals (q)
        % q1 = (100 - 90)/max(90, 10) = 10/90 = 0.1111
        % q2 = (80 - 80)/max(80, 10) = 0
        % q3 = (200 - 220)/max(220, 10) = -20/220 = -0.0909
        % all_q = [0.1111, 0, -0.0909]
        % abs(all_q) = [0.1111, 0, 0.0909]
        % median(abs(all_q)) = 0.0909...
        % EXPECTED gamma: 1.4826 * 0.0909 = 0.1348
        
        expectedGamma = 1.4826 * median(abs([0.11111111; 0; -0.09090909]));
        testCase.verifyEqual(model.gamma, expectedGamma, 'AbsTol', 1e-4);
    end


    function testEmptyInput(testCase)
        % Test empty input exceptions
        testCase.verifyError(@() CMS2QuantSolver.estimateDatasetNoiseModel([]), ...
            'CLogger:LoggedError');
    end

end
end
