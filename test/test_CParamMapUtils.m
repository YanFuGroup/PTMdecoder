function tests = test_CParamMapUtils
% TEST_CPARAMMAPUTILS Unit tests for CParamMapUtils helper methods.
% Input:
%   (none)
% Output:
%   tests (matlab.unittest.Test)

tests = functiontests(localfunctions);
end

function testParseQuotedListBasicSemicolon(testCase)
% TESTPARSEQUOTEDLISTBASICSEMICOLON Parse quoted list with ';' delimiter.

raw = '"alpha";"beta";"gamma"';
actual = CParamMapUtils.parseQuotedList(raw, ';');
expected = {'alpha', 'beta', 'gamma'};
testCase.verifyEqual(actual, expected);
end

function testParseQuotedListDefaultDelimiter(testCase)
% TESTPARSEQUOTEDLISTDEFAULTDELIMITER Parse quoted list using default delimiter.

raw = '"left";"right"';
actual = CParamMapUtils.parseQuotedList(raw);
expected = {'left', 'right'};
testCase.verifyEqual(actual, expected);
end

function testParseQuotedListCustomDelimiter(testCase)
% TESTPARSEQUOTEDLISTCUSTOMDELIMITER Parse quoted list with custom delimiter.

raw = '"x"|"y"|"z"';
actual = CParamMapUtils.parseQuotedList(raw, '|');
expected = {'x', 'y', 'z'};
testCase.verifyEqual(actual, expected);
end

function testParseQuotedListIgnoresUnquotedSegments(testCase)
% TESTPARSEQUOTEDLISTIGNORESUNQUOTEDSEGMENTS Keep only quoted segments.

raw = 'abc;"keep1";;"keep2";123';
actual = CParamMapUtils.parseQuotedList(raw, ';');
expected = {'keep1', 'keep2'};
testCase.verifyEqual(actual, expected);
end

function testParseQuotedListSupportsStringScalar(testCase)
% TESTPARSEQUOTEDLISTSUPPORTSSTRINGSCALAR Support MATLAB string scalar input.

raw = """s1"";""s2""";
actual = CParamMapUtils.parseQuotedList(raw, ';');
expected = {'s1', 's2'};
testCase.verifyEqual(actual, expected);
end

function testParseQuotedListQuotedSegmentsWithSpaces(testCase)
% TESTPARSEQUOTEDLISTQUOTEDSEGMENTSWITHSPACES Parse quoted segments with spaces around delimiters.

raw = ' "a" ; "b" ';
actual = CParamMapUtils.parseQuotedList(raw, ';');
expected = {'a', 'b'};
testCase.verifyEqual(actual, expected);
end

function testParseQuotedListEmptyInput(testCase)
% TESTPARSEQUOTEDLISTEMPTYINPUT Empty input returns empty cell array.

actual = CParamMapUtils.parseQuotedList('', ';');
expected = {};
testCase.verifyEqual(actual, expected);
end

function testParseQuotedListInvalidTypeError(testCase)
% TESTPARSEQUOTEDLISTINVALIDTYPEERROR Invalid raw value type should error.

testCase.verifyError(@() CParamMapUtils.parseQuotedList(123, ';'), ...
    'CParamMapUtils:InvalidQuotedListType');
end

function testGetOptionalLogicalParsesCommonValues(testCase)
% TESTGETOPTIONALLOGICALPARSESCOMMONVALUES Parse logical-like values from map.

param_map = containers.Map('KeyType', 'char', 'ValueType', 'any');
param_map('flag_true') = 'true';
param_map('flag_false') = '0';

actual_true = CParamMapUtils.getOptionalLogical(param_map, 'flag_true', false);
actual_false = CParamMapUtils.getOptionalLogical(param_map, 'flag_false', true);
default_value = CParamMapUtils.getOptionalLogical(param_map, 'missing_key', true);

testCase.verifyTrue(actual_true);
testCase.verifyFalse(actual_false);
testCase.verifyTrue(default_value);
end

function testGetOptionalLogicalInvalidValueError(testCase)
% TESTGETOPTIONALLOGICALINVALIDVALUEERROR Invalid logical-like value should error.

param_map = containers.Map('KeyType', 'char', 'ValueType', 'any');
param_map('bad_flag') = 'not_a_boolean';

testCase.verifyError(@() CParamMapUtils.getOptionalLogical(param_map, 'bad_flag', false), ...
    'CParamMapUtils:InvalidLogicalParamValue');
end
