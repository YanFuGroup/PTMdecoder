classdef test_CPrecursorProfileUtils < matlab.unittest.TestCase
    % Unit tests for spectrum-name parsing utilities.

    methods (Test)
        function testParseMS2ScanNumberFromDottedName(testCase)
            % Parse second token from dotted spectrum name.
            % Inputs:
            %    testCase (matlab.unittest.TestCase)
            %        Test case context.
            % Outputs:
            %    (none)

            scan_number = CMS2SpecNameUtils.parseMS2ScanNumber('scan.12345.12345.2');
            testCase.verifyEqual(scan_number, 12345);
        end


        function testParseMS2ScanNumberFromRawNumber(testCase)
            % Parse raw scan number when spectrum name has no dot.
            % Inputs:
            %    testCase (matlab.unittest.TestCase)
            %        Test case context.
            % Outputs:
            %    (none)

            scan_number = CMS2SpecNameUtils.parseMS2ScanNumber('67890');
            testCase.verifyEqual(scan_number, 67890);
        end


        function testParseMS2ScanNumberRejectsInvalidName(testCase)
            % Throw a clear error for non-numeric spectrum name.
            % Inputs:
            %    testCase (matlab.unittest.TestCase)
            %        Test case context.
            % Outputs:
            %    (none)

            testCase.verifyError(@() CMS2SpecNameUtils.parseMS2ScanNumber('scan.onlytext'), ...
                'CMS2SpecNameUtils:InvalidSpectrumName');
        end
    end
end
