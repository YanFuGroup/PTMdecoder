classdef test_CMS2PeakMatcher < matlab.unittest.TestCase

    methods (Test)
        function testMatchAggregatesAndFilters(testCase)
            expPeaks = [
                100.05, 10;
                100.08, 5;
                101.00, 20;
                150.00, 99
            ];
            vNonRedunTheoryIonMz = [
                100.00;
                101.00;
                200.00
            ];

            matched = CMS2PeakMatcher.match(expPeaks, vNonRedunTheoryIonMz, 0.1);

            testCase.verifyEqual(matched, [1, 15; 2, 20], 'AbsTol', 1e-12);
        end

        function testProcessPeaksNormalizeAndThreshold(testCase)
            expPeaks = [
                10, 1;
                20, 5;
                30, 10
            ];

            processed = CMS2PeakMatcher.processPeaks(expPeaks, 0.2);

            expected = [
                20, 0.5, 5;
                30, 1.0, 10
            ];
            testCase.verifyEqual(processed, expected, 'AbsTol', 1e-12);
        end

        function testPostMatchFilterPeptidoformsRemovesStrictSubsets(testCase)
            matchedExpPeaks = [
                1, 0.9, 90;
                2, 0.8, 80;
                3, 0.7, 70
            ];

            massArrangement = [
                10, 20;
                30, 40;
                50, 60
            ];

            baseInfo = [
                500, 1, 1, 1, 0, 1;
                600, 1, 2, 1, 0, 2;
                700, 2, 3, 1, 0, 3
            ];

            imp1 = [1; 1; 0];
            imp2 = [1; 1; 1];
            imp3 = [0; 1; 0];
            vNonRedunTheoryIonMz = [baseInfo, imp1, imp2, imp3];

            [massOut, ionsOut] = CMS2PeakMatcher.postMatchFilterPeptidoforms( ...
                matchedExpPeaks, massArrangement, vNonRedunTheoryIonMz);

            testCase.verifyEqual(massOut, [30, 40], 'AbsTol', 1e-12);
            testCase.verifyEqual(size(ionsOut), [3, 7]);
            testCase.verifyEqual(ionsOut(:, 7), imp2, 'AbsTol', 1e-12);
        end
    end
end
