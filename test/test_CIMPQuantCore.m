function tests = test_CIMPQuantCore
% Validate CIMPQuantCore threshold validation and early diagnostics.
tests = functiontests(localfunctions);
end

function testQuantGroupRejectsInvalidMinimumBeforeDatasetAccess(testCase)
invalid_values = {0, -1, 1.5, NaN, Inf};
for idx = 1:numel(invalid_values)
    testCase.verifyError(@() CIMPQuantCore.quantGroup( ...
        [], 'raw.mgf', 1, 1, 1, 500, 501, 2, 1, 0.01, 0, invalid_values{idx}), ...
        'CIMPQuantCore:InvalidMinXicNonzeroPoints');
end
end

function testQuantGroupDefaultsAndReportsInsufficientInputs(testCase)
[has_nonzero_imp, ~, ~, ~, ~, ~, diagnostics] = CIMPQuantCore.quantGroup( ...
    [], 'raw.mgf', 1, 1, 1, 500, 501, 2, 2, 0.01, 0);

testCase.verifyFalse(has_nonzero_imp);
testCase.verifyEqual(diagnostics.reason, 'insufficient_psm_inputs');
testCase.verifyEqual(diagnostics.min_nonzero_points, 5);
testCase.verifyEqual(diagnostics.candidate_peak_count, 0);
testCase.verifyEqual(diagnostics.filtered_sparse_peak_count, 0);
end
