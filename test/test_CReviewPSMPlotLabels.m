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

function testPlotMatchColorsSingleImpPeakWhenConfigured(testCase)
review = makeReviewPsm(makePeptides({'IMP06'}));
review.m_all_match_ions = [1, 5, 1, 1, 1];
options = makeColorOptions({'IMP06'}, {[0.10, 0.20, 0.30]}, [1]);

review.plotMatch(options);

colored_stems = findall(gca, 'Tag', 'CReviewPSMImpColoredStem');
testCase.assertNumElements(colored_stems, 1);
testCase.verifyEqual(colored_stems.Color, [0.10, 0.20, 0.30], 'AbsTol', 1e-12);
testCase.verifyEqual(colored_stems.YData, [0, 100], 'AbsTol', 1e-12);
end

function testPlotMatchStacksSharedPeakByUniqueImpProportions(testCase)
review = makeReviewPsm(makePeptides({'IMP06', 'IMP08'}));
review.m_all_match_ions = [
    1, 5, 1, 1, 1;
    2, 7, 1, 1, 1;
    1, 5, 1, 1, 2
];
options = makeColorOptions( ...
    {'IMP06', 'IMP08'}, ...
    {[0.10, 0.20, 0.30], [0.70, 0.80, 0.90]}, ...
    [3, 1]);

review.plotMatch(options);

segments = findall(gca, 'Tag', 'CReviewPSMImpColoredSegment');
testCase.assertNumElements(segments, 2);
segment_colors = vertcat(segments.Color);
testCase.verifyTrue(any(all(abs(segment_colors - [0.10, 0.20, 0.30]) < 1e-12, 2)));
testCase.verifyTrue(any(all(abs(segment_colors - [0.70, 0.80, 0.90]) < 1e-12, 2)));
segment_y = sortrows(vertcat(segments.YData), 1);
testCase.verifyEqual(segment_y, [0, 75; 75, 100], 'AbsTol', 1e-12);
end

function testPlotMatchFallsBackWhenColorConfigurationIsIncomplete(testCase)
review = makeReviewPsm(makePeptides({'IMP06', 'IMP08'}));
review.m_all_match_ions = [
    1, 5, 1, 1, 1;
    1, 5, 1, 1, 2
];
options = makeColorOptions({'IMP06'}, {[0.10, 0.20, 0.30]}, [1]);

review.plotMatch(options);

colored_segments = findall(gca, 'Tag', 'CReviewPSMImpColoredSegment');
colored_stems = findall(gca, 'Tag', 'CReviewPSMImpColoredStem');
testCase.verifyEmpty(colored_segments);
testCase.verifyEmpty(colored_stems);

labels = getFigureTextLabels();
testCase.verifyTrue(any(strcmp(labels, 'IMP06,IMP08:y5+')));
end

function testPlotMatchRejectsInvalidColorOptions(testCase)
review = makeReviewPsm(makePeptides({'IMP06'}));
review.m_all_match_ions = [1, 5, 1, 1, 1];
options = makeColorOptions({'IMP06'}, {[1.20, 0.20, 0.30]}, [1]);

testCase.verifyError(@() review.plotMatch(options), 'CReviewPSM:InvalidImpColor');
end

function testPlotMatchRejectsInvalidProportionOptions(testCase)
review = makeReviewPsm(makePeptides({'IMP06'}));
review.m_all_match_ions = [1, 5, 1, 1, 1];
options = makeColorOptions({'IMP06'}, {[0.10, 0.20, 0.30]}, [NaN]);

testCase.verifyError(@() review.plotMatch(options), 'CReviewPSM:InvalidImpProportion');
end

function testPlotMatchFallsBackForPeakWithOnlyZeroProportionImps(testCase)
review = makeReviewPsm(makePeptides({'IMP06', 'IMP10'}));
review.m_all_match_ions = [
    1, 5, 1, 1, 1;
    1, 6, 1, 2, 2
];
options = makeColorOptions( ...
    {'IMP06', 'IMP10'}, ...
    {[0.10, 0.20, 0.30], [0.70, 0.80, 0.90]}, ...
    [1, 0]);

review.plotMatch(options);

colored_stems = findall(gca, 'Tag', 'CReviewPSMImpColoredStem');
testCase.verifyNumElements(colored_stems, 1);
legacy_zero_sum_peak = findall(gca, 'Tag', 'CReviewPSMLegacyZeroProportionPeak');
testCase.assertNumElements(legacy_zero_sum_peak, 1);
testCase.verifyEqual(legacy_zero_sum_peak.Color, [1, 0, 0], 'AbsTol', 1e-12);
testCase.verifyEqual(legacy_zero_sum_peak.XData, [200, 200], 'AbsTol', 1e-12);
testCase.verifyEqual(legacy_zero_sum_peak.YData, [0, 80], 'AbsTol', 1e-12);
end

function testPlotMatchShowsOutsideImpColorLegendWhenRequested(testCase)
review = makeReviewPsm(makePeptides({'IMP06', 'IMP08'}));
review.m_all_match_ions = [
    1, 5, 1, 1, 1;
    1, 5, 1, 2, 2
];
options = makeColorOptions( ...
    {'IMP06', 'IMP08', 'IMP10'}, ...
    {[0.10, 0.20, 0.30], [0.70, 0.80, 0.90], [0.30, 0.40, 0.50]}, ...
    [3, 1, 2]);
options.show_imp_color_legend = true;

review.plotMatch(options);

legends = findall(gcf, 'Type', 'legend', 'Tag', 'CReviewPSMImpColorLegend');
testCase.assertNumElements(legends, 1);
testCase.verifyEqual(string(legends.String), string({'IMP06', 'IMP08'}));
testCase.verifyEqual(char(legends.Location), 'eastoutside');
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

function options = makeColorOptions(labels, colors, proportions)
color_map = containers.Map();
proportion_map = containers.Map();
for idx = 1:numel(labels)
    color_map(labels{idx}) = colors{idx};
    proportion_map(labels{idx}) = proportions(idx);
end
options = struct( ...
    'imp_colors', color_map, ...
    'imp_proportions', proportion_map);
end

function labels = getFigureTextLabels()
text_objects = findall(gcf, 'Type', 'text');
labels = get(text_objects, 'String');
if ischar(labels) || isstring(labels)
    labels = {char(labels)};
end
labels = cellfun(@char, labels, 'UniformOutput', false);
end
