function tests = test_CReviewPSMPlotLabels
% TEST_CREVIEWPSMPLOTLABELS Validate CReviewPSM peak label formatting.

tests = functiontests(localfunctions);
end

function teardown(~)
close all force;
end

function testPlotMatchFallsBackToLocalPeptideIndex(testCase)
review = makeReviewPsm(makePeptides({'', ''}));
review.m_all_match_ions = [1, 5, 1, 1, 2];

review.plotMatch();

labels = getFigureTextLabels();
testCase.verifyTrue(any(strcmp(labels, '2:y5+')));
end

function testPlotMatchUsesPeptideLabelWhenAvailable(testCase)
review = makeReviewPsm(makePeptides({'IMP06', 'IMP07'}));
review.m_all_match_ions = [2, 4, 1, 1, 1];

review.plotMatch();

labels = getFigureTextLabels();
testCase.verifyTrue(any(strcmp(labels, 'IMP06:b4+')));
end

function testPlotMatchDefaultShowsTwoMatchesAndEllipsis(testCase)
review = makeReviewPsm(makePeptides({'IMP06', 'IMP07', 'IMP08', 'IMP10'}));
review.m_all_match_ions = [
    1, 5, 1, 1, 1;
    1, 5, 1, 1, 2;
    2, 7, 1, 1, 3;
    1, 8, 1, 1, 4
];

review.plotMatch();

labels = getFigureTextLabels();
testCase.verifyTrue(any(strcmp(labels, 'IMP06,IMP07:y5+;IMP08:b7+;...')));
end

function testPlotMatchCanShowFourMatches(testCase)
review = makeReviewPsm(makePeptides({'IMP06', 'IMP07', 'IMP08', 'IMP10'}));
review.m_all_match_ions = [
    1, 5, 1, 1, 1;
    1, 5, 1, 1, 2;
    2, 7, 1, 1, 3;
    2, 8, 1, 1, 4
];

review.plotMatch(struct('max_labels_per_peak', 4));

labels = getFigureTextLabels();
testCase.verifyTrue(any(strcmp(labels, 'IMP06,IMP07:y5+;IMP08:b7+;IMP10:b8+')));
end

function testPlotMatchGroupsFallbackPeptideIndexes(testCase)
review = makeReviewPsm(makePeptides({'', '', ''}));
review.m_all_match_ions = [
    1, 5, 1, 1, 1;
    1, 5, 1, 1, 2;
    1, 5, 1, 1, 3
];

review.plotMatch(struct('max_labels_per_peak', 4));

labels = getFigureTextLabels();
testCase.verifyTrue(any(strcmp(labels, '1,2,3:y5+')));
end

function testPlotMatchUsesPercentScaleAndTopMargin(testCase)
review = makeReviewPsm(makePeptides({'IMP06'}));
review.m_all_match_ions = [1, 5, 1, 1, 1];

review.plotMatch();

ax = gca;
ylim_values = ylim(ax);
testCase.verifyGreaterThan(max(ylim_values), 100);

testCase.verifyEqual(ax.YAxis.Visible, matlab.lang.OnOffSwitchState.off);

custom_y_axis = findall(ax, 'Tag', 'CReviewPSMCustomYAxisLine');
testCase.verifyNumElements(custom_y_axis, 1);
testCase.verifyEqual(custom_y_axis.YData, [0, 100], 'AbsTol', 1e-12);

custom_tick_labels = findall(ax, 'Tag', 'CReviewPSMCustomYAxisTickLabel');
tick_label_strings = sort(string(get(custom_tick_labels, 'String')));
testCase.verifyEqual(tick_label_strings, string({'0'; '100'; '20'; '40'; '60'; '80'}));
end

function testPlotMatchShowsPrecursorByDefault(testCase)
review = makeReviewPsm(makePeptides({'IMP06'}));
review.m_all_match_ions = [1, 5, 1, 1, 1];

review.plotMatch();

labels = getFigureTextLabels();
testCase.verifyTrue(any(strcmp(labels, '[M]++')));
end

function testPlotMatchCanHidePrecursorLabel(testCase)
review = makeReviewPsm(makePeptides({'IMP06'}));
review.m_all_match_ions = [1, 5, 1, 1, 1];

review.plotMatch(struct('show_precursor', false));

labels = getFigureTextLabels();
testCase.verifyFalse(any(strcmp(labels, '[M]++')));
end

function review = makeReviewPsm(peptides)
spectrum = struct( ...
    'peaks', [100, 100; 200, 80; 300, 70], ...
    'pre_charge', 2, ...
    'pre_mz', 500 ...
);
tolerance = struct('value', 0.02, 'is_ppm', false);
review = CReviewPSM(peptides, spectrum, tolerance);
end

function peptides = makePeptides(labels)
base = struct('seq', 'AA', 'mod_mass', [], 'mod_pos', []);
peptides = repmat(base, 1, numel(labels));
for idx = 1:numel(labels)
    if ~isempty(labels{idx})
        peptides(idx).label = labels{idx};
    end
end
end

function labels = getFigureTextLabels()
text_objects = findall(gcf, 'Type', 'text');
labels = get(text_objects, 'String');
if ischar(labels) || isstring(labels)
    labels = {char(labels)};
end
labels = cellfun(@char, labels, 'UniformOutput', false);
end
