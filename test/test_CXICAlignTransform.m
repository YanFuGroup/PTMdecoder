classdef test_CXICAlignTransform < matlab.unittest.TestCase
    % Unit tests for CXICAlignTransform

    methods (Test)
        function testFitLinearModel(testCase)
            ref_rts = (1:10)';
            target_rts = 2 .* ref_rts + 3;
            options = struct('outlier_k', 0, 'num_bins', 1);

            model = CXICAlignTransform.fit(ref_rts, target_rts, options);

            testCase.verifyTrue(model.has_model);
            testCase.verifyEqual(model.slope, 2, 'AbsTol', 1e-10);
            testCase.verifyEqual(model.intercept, 3, 'AbsTol', 1e-10);
        end

        function testFitOutlierRemovalMad(testCase)
            ref_rts = (1:10)';
            target_rts = 2 .* ref_rts + 3;
            target_rts(end) = 100; % outlier

            model_no = CXICAlignTransform.fit(ref_rts, target_rts, struct('outlier_k', 0));
            model_yes = CXICAlignTransform.fit(ref_rts, target_rts, struct('outlier_k', 2, 'outlier_method', 'mad'));

            dist_no = abs(model_no.slope - 2);
            dist_yes = abs(model_yes.slope - 2);
            testCase.verifyLessThan(dist_yes, dist_no);
            testCase.verifyNotEmpty(model_yes.rt_residual_threshold);
        end

        function testFitInsufficientPoints(testCase)
            ref_rts = 1;
            target_rts = 2;
            model = CXICAlignTransform.fit(ref_rts, target_rts, struct());

            testCase.verifyFalse(model.has_model);
            testCase.verifyEqual(model.slope, 1, 'AbsTol', 1e-12);
            testCase.verifyEqual(model.intercept, 0, 'AbsTol', 1e-12);
        end

        function testFitBins(testCase)
            ref_rts = (0:8)';
            target_rts = ref_rts + [0; 0; 0; 2; 2; 2; -1; -1; -1];
            options = struct('num_bins', 3, 'min_per_bin', 3, 'outlier_k', 0);

            model = CXICAlignTransform.fit(ref_rts, target_rts, options);

            testCase.verifyTrue(model.has_model);
            testCase.verifyEqual(numel(model.bin_centers), 3);
            testCase.verifyEqual(numel(model.bin_offsets), 3);
            testCase.verifyTrue(any(abs(model.bin_offsets) > 0.05));
        end

        function testApplyLinearNoBins(testCase)
            model = struct('slope', 2, 'intercept', 1, 'bin_centers', [], ...
                'bin_offsets', [], 'bin_min', [], 'bin_max', []);
            rt_in = [1; 2; 3];

            rt_out = CXICAlignTransform.apply(model, rt_in);

            testCase.verifyEqual(rt_out, [3; 5; 7], 'AbsTol', 1e-12);
        end

        function testApplyWithBinsClamp(testCase)
            model = struct('slope', 1, 'intercept', 0, ...
                'bin_centers', [2, 4], 'bin_offsets', [1, -1], ...
                'bin_min', 2, 'bin_max', 4);
            rt_in = [1; 2; 3; 5];

            rt_out = CXICAlignTransform.apply(model, rt_in);

            testCase.verifyEqual(rt_out, [2; 3; 3; 4], 'AbsTol', 1e-12);
        end
    end
end
