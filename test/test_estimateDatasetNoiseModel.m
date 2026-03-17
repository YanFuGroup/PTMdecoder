classdef test_estimateDatasetNoiseModel < matlab.unittest.TestCase
% Unit test for CMS2QuantSolver.estimateDatasetNoiseModel

methods (Test)

    function testNormalDataset(testCase)
        % Test basic pooled MLE and MAD calculation behavior
        
        % Setup mock inputs
        inputs(1).filteredOutExpPeakCount = 100;
        inputs(1).filteredOutExpPeakSqSum = 0.04; % pooled variance = 4e-4 => sigma = 0.02
        inputs(1).matchedExpPeaks = [1, 0.90, 900; 2, 0.80, 800]; % [ionIndex, normInt, compatibility]
        inputs(1).fittedMatchedPeakIntensities = [0.81; 0.80]; % First has residual 0.09, second 0
        
        inputs(2).filteredOutExpPeakCount = 50;
        inputs(2).filteredOutExpPeakSqSum = 0.02; % Total count=150, SqSum=0.06 => sigma=0.02
        inputs(2).matchedExpPeaks = [1, 0.95, 950; 2, 0.05, 50];
        inputs(2).fittedMatchedPeakIntensities = [0.99; 0.05]; % First has residual -0.04, second 0
        
        % Execute
        model = CMS2QuantSolver.estimateDatasetNoiseModel(inputs);
        
        % Verify sigma_base
        testCase.verifyEqual(model.sigma_base, 0.02, 'AbsTol', 1e-6);
        testCase.verifyEqual(model.tau_floor, 0.1, 'AbsTol', 1e-6);
        
        % Theoretical high signals: 0.81, 0.80, 0.99 (all > 0.1)
        % Normalized y: 0.90, 0.80, 0.95
        % Residuals (q)
        % q1 = (0.90 - 0.81)/max(0.81, 0.1) = 0.09/0.81 = 0.1111
        % q2 = (0.80 - 0.80)/max(0.80, 0.1) = 0
        % q3 = (0.95 - 0.99)/max(0.99, 0.1) = -0.04/0.99 = -0.0404
        % all_q = [0.1111, 0, -0.0404]
        % abs(all_q) = [0.1111, 0, 0.0404]
        % median(abs(all_q)) = 0.0404...
        % EXPECTED gamma: 1.4826 * 0.0404 = 0.0599
        
        expectedGamma = 1.4826 * median(abs([0.11111111; 0; -0.04040404]));
        testCase.verifyEqual(model.gamma, expectedGamma, 'AbsTol', 1e-6);
    end


    function testEmptyInput(testCase)
        % Test empty input exceptions
        testCase.verifyError(@() CMS2QuantSolver.estimateDatasetNoiseModel([]), ...
            'CLogger:LoggedError');
    end

end
end
