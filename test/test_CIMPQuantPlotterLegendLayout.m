function tests = test_CIMPQuantPlotterLegendLayout
% Validate adaptive right-side legend layout for XIC plots.
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
test_dir = fileparts(mfilename('fullpath'));
project_dir = fileparts(test_dir);
addpath(project_dir);
addpath(fullfile(test_dir, 'helpers'));
testCase.TestData.output_dir = fullfile(test_dir, 'output', 'legend_layout');
if ~exist(testCase.TestData.output_dir, 'dir')
    mkdir(testCase.TestData.output_dir);
end
end

function teardown(testCase)
close all force;
end

function testEmptyXicLayoutOverrideReturnsDefaultLayout(testCase)
layout = CIMPQuantPlotter.resolveXicLegendLayoutConfig([]);
default_layout = CIMPQuantPlotter.getXicLegendLayoutConfig();

testCase.verifyEqual(layout, default_layout);
end

function testPartialXicLayoutOverrideMergesDefaults(testCase)
override = struct();
override.figure_width_px = 1333;
override.figure_height_px = 667;
override.axes_width_fraction = 0.75;

layout = CIMPQuantPlotter.resolveXicLegendLayoutConfig(override);
default_layout = CIMPQuantPlotter.getXicLegendLayoutConfig();

testCase.verifyEqual(layout.figure_width_px, 1333);
testCase.verifyEqual(layout.figure_height_px, 667);
testCase.verifyEqual(layout.axes_width_fraction, 0.75);
testCase.verifyEqual(layout.left_margin_px, default_layout.left_margin_px);
testCase.verifyTrue(isfield(layout, 'axes_width_px'));
testCase.verifyTrue(isfield(layout, 'legend_max_width_px'));

available_width_px = layout.figure_width_px ...
    - layout.left_margin_px - layout.axes_legend_gap_px ...
    - layout.right_margin_px - layout.legend_print_safety_px;
testCase.verifyEqual(layout.axes_width_px, round(available_width_px * layout.axes_width_fraction));
testCase.verifyEqual(layout.legend_max_width_px, available_width_px - layout.axes_width_px);
end

function testInvalidXicLayoutOverrideFailsClearly(testCase)
testCase.verifyError(@() CIMPQuantPlotter.resolveXicLegendLayoutConfig('figure_width_px=1333'), ...
    'CIMPQuantPlotter:InvalidXicLayoutOverride');
end

function testLegendLabelsDoNotUseXicOfOrTexSubscripts(testCase)
layout = CIMPQuantPlotter.getXicLegendLayoutConfig();
labels = {'_PEPT{phospho}IDESK_', '_PEPTIDES{phospho}K_'};
wrapped = CIMPQuantPlotter.prepareXicLegendLabels(labels, layout);

testCase.verifyFalse(any(contains(wrapped, 'XIC of')));
testCase.verifyEqual(wrapped{1}, '_PEPT{phospho}IDESK_');
testCase.verifyEqual(wrapped{2}, '_PEPTIDES{phospho}K_');
end

function testLongLegendWrapsWithinConfiguredLimit(testCase)
layout = CIMPQuantPlotter.getXicLegendLayoutConfig();
layout.legend_max_width_px = 210;
labels = {'_PEPT{phospho}IDESK_WITH_EXTRA_LONG_CONTEXT_FOR_LAYOUT_'};
wrapped = CIMPQuantPlotter.prepareXicLegendLabels(labels, layout);

testCase.verifyTrue(contains(wrapped{1}, newline));
testCase.verifyFalse(contains(wrapped{1}, 'XIC of'));
testCase.verifyTrue(contains(wrapped{1}, '{phospho}'));
end

function testLegendWrapsByFixedCharacterCount(testCase)
layout = CIMPQuantPlotter.getXicLegendLayoutConfig();
layout.legend_max_width_px = 130;
labels = {'_TRTDSTS{Phospho}DGRPAWMR_'};
wrapped = CIMPQuantPlotter.prepareXicLegendLabels(labels, layout);

testCase.verifyEqual(wrapped{1}, sprintf('_TRTDSTS{Pho\nspho}DGRPAWM\nR_'));
end

function testLegendWidthControlsAutomaticLineLength(testCase)
layout = CIMPQuantPlotter.getXicLegendLayoutConfig();
label = '_TRTDSTS{Phospho}DGRPAWMR_WITH_EXTRA_SUFFIX_';

layout.legend_max_width_px = 180;
wrapped_narrow = CIMPQuantPlotter.prepareXicLegendLabels({label}, layout);

layout.legend_max_width_px = 360;
wrapped_wide = CIMPQuantPlotter.prepareXicLegendLabels({label}, layout);

testCase.verifyGreaterThan(numel(splitlines(wrapped_narrow{1})), numel(splitlines(wrapped_wide{1})));
testCase.verifyEqual(erase(wrapped_narrow{1}, newline), label);
testCase.verifyEqual(erase(wrapped_wide{1}, newline), label);
end

function testLongLegendDiagnosticsPreventClipping(testCase)
layout = CIMPQuantPlotter.getXicLegendLayoutConfig();
layout.figure_width_px = 1057;
layout.axes_width_fraction = 4 / 7;
layout.figure_height_px = 460;
layout.left_margin_px = 115;
layout.bottom_margin_px = 105;
layout.top_margin_px = 50;
layout.axes_legend_gap_px = 40;
layout.right_margin_px = 30;
layout.legend_min_width_px = 160;
layout = finalizeXicLegendTestLayout(layout);
testCase.verifyEqual(layout.figure_width_px, 1057);
testCase.verifyEqual(layout.axes_width_px, 480);
testCase.verifyEqual(layout.legend_max_width_px, 360);

output_base = fullfile(testCase.TestData.output_dir, 'long_legend_export');
for ext = {'.png', '.pdf'}
    out_file = [output_base, ext{1}];
    if exist(out_file, 'file')
        delete(out_file);
    end
end

rt = linspace(10, 12, 80);
total_xic = {rt, exp(-((rt - 11).^2) / 0.08)};
ric = {
    rt, 0.60 .* total_xic{2};
    rt, 0.40 .* total_xic{2}
};
names = {
    '_PEPT{phospho}IDESK_WITH_EXTRA_LONG_CONTEXT_FOR_LAYOUT_';
    '_PEPTIDES{phospho}K_WITH_EXTRA_LONG_CONTEXT_FOR_LAYOUT_'
};
colors = [
    0.0000, 0.4470, 0.7410;
    0.8500, 0.3250, 0.0980
];

[color_map, legend_map] = buildXicLegendMaps(names, colors);
diagnostics = CIMPQuantPlotter.plotXicGroupWithLayout( ...
    ric, total_xic, [10, 12], names, output_base, color_map, legend_map, layout);

testCase.verifyTrue(diagnostics.is_legend_inside_figure, diagnostics.summary);
testCase.verifyTrue(diagnostics.is_text_inside_figure, diagnostics.summary);
testCase.verifyTrue(diagnostics.is_text_clear_of_axes, diagnostics.summary);
testCase.verifyGreaterThan(diagnostics.min_text_right_margin_px, 0, diagnostics.summary);

fig_w = diagnostics.figure_position_px(3);
fig_h = diagnostics.figure_position_px(4);
testCase.verifyGreaterThanOrEqual(fig_w, 700);
testCase.verifyLessThanOrEqual(fig_w, 1200);
testCase.verifyGreaterThanOrEqual(fig_w / fig_h, 2);
testCase.verifyLessThanOrEqual(fig_w / fig_h, 4);

text_extents = diagnostics.legend_text_extents_px;
testCase.verifyGreaterThan(size(text_extents, 1), 0, 'Legend text extents are empty.');
for idx_extent = 1:size(text_extents, 1)
    text_left = text_extents(idx_extent, 1);
    text_bottom = text_extents(idx_extent, 2);
    text_right = text_extents(idx_extent, 1) + text_extents(idx_extent, 3);
    text_top = text_extents(idx_extent, 2) + text_extents(idx_extent, 4);
    testCase.verifyGreaterThanOrEqual(text_left, 0, diagnostics.summary);
    testCase.verifyGreaterThanOrEqual(text_bottom, 0, diagnostics.summary);
    testCase.verifyLessThanOrEqual(text_right, fig_w, diagnostics.summary);
    testCase.verifyLessThanOrEqual(text_top, fig_h, diagnostics.summary);
end

plot_band_center_y = diagnostics.axes_center_y_px;
text_center_y = diagnostics.legend_center_y_px;
vertical_center_tolerance_px = 2;
testCase.verifyLessThan(abs(text_center_y - plot_band_center_y), vertical_center_tolerance_px, ...
    sprintf('Legend text block not vertically centered: text_center=%.1f, band_center=%.1f', ...
    text_center_y, plot_band_center_y));

metrics = assertXicLegendExports(output_base, layout, diagnostics);
testCase.verifyGreaterThan(metrics.legend_right_margin_px, 0);
testCase.verifyLessThanOrEqual(metrics.right_blank_fraction, 0.20);
end

function testDelimiterFreeLongLegendUsesHardWrapAndRenderedProbe(testCase)
layout = CIMPQuantPlotter.getXicLegendLayoutConfig();
layout.figure_width_px = 1007;
layout.axes_width_fraction = 0.65;
layout.figure_height_px = 480;
layout.left_margin_px = 115;
layout.bottom_margin_px = 110;
layout.top_margin_px = 55;
layout.axes_legend_gap_px = 36;
layout.right_margin_px = 24;
layout.legend_min_width_px = 150;
layout = finalizeXicLegendTestLayout(layout);
testCase.verifyEqual(layout.axes_width_px, 520);
testCase.verifyEqual(layout.legend_max_width_px, 280);

output_base = fullfile(testCase.TestData.output_dir, 'delimiter_free_hard_wrap');
for ext = {'.png', '.pdf', '.svg'}
    out_file = [output_base, ext{1}];
    if exist(out_file, 'file')
        delete(out_file);
    end
end

rt = linspace(10, 12, 90);
total_xic = {rt, exp(-((rt - 11).^2) / 0.08)};
ric = {
    rt, 0.60 .* total_xic{2};
    rt, 0.40 .* total_xic{2}
};
names = {
    'PEPTPHOSPHOIDESKAAAAAAAAAAAAAAAAAAAAAAAAAAAA';
    'PEPTIDESPHOSPHOKBBBBBBBBBBBBBBBBBBBBBBBBBBBB'
};
colors = [
    0.0000, 0.4470, 0.7410;
    0.8500, 0.3250, 0.0980
];

[color_map, legend_map] = buildXicLegendMaps(names, colors);
diagnostics = CIMPQuantPlotter.plotXicGroupWithLayout( ...
    ric, total_xic, [10, 12], names, output_base, color_map, legend_map, layout);

testCase.verifyGreaterThanOrEqual(diagnostics.rendered_legend_right_margin_px, 8, diagnostics.summary);
testCase.verifyGreaterThan(diagnostics.final_legend_line_chars, 0);
testCase.verifyTrue(any(contains(diagnostics.final_legend_labels, newline)), ...
    'Delimiter-free long labels should be hard-wrapped.');
compact_labels = cell(size(diagnostics.final_legend_labels));
for idx_label = 1:numel(diagnostics.final_legend_labels)
    compact_labels{idx_label} = erase(diagnostics.final_legend_labels{idx_label}, newline);
end
testCase.verifyEqual(sort(compact_labels(:)), sort(names(:)));

metrics = assertXicLegendExports(output_base, layout, diagnostics);
testCase.verifyGreaterThanOrEqual(metrics.legend_right_margin_px, 8);
testCase.verifyTrue(metrics.svg_checked);
end

function testUnsatisfiableFixedWidthLegendFailsClearly(testCase)
layout = CIMPQuantPlotter.getXicLegendLayoutConfig();
layout.figure_width_px = 674;
layout.axes_width_fraction = 480 / (480 + 28);
layout.figure_height_px = 260;
layout.left_margin_px = 80;
layout.bottom_margin_px = 50;
layout.top_margin_px = 40;
layout.axes_legend_gap_px = 30;
layout.right_margin_px = 24;
layout.legend_min_width_px = 28;
layout.legend_padding_px = 16;
layout = finalizeXicLegendTestLayout(layout);
testCase.verifyEqual(layout.axes_width_px, 480);
testCase.verifyEqual(layout.legend_max_width_px, 28);

output_base = fullfile(testCase.TestData.output_dir, 'unsatisfiable_fixed_width');
writeStaleXicArtifacts(testCase, output_base);
rt = linspace(10, 12, 40);
total_xic = {rt, exp(-((rt - 11).^2) / 0.08)};
ric = {rt, total_xic{2}};
names = {'WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW'};
colors = [0.0000, 0.4470, 0.7410];

[color_map, legend_map] = buildXicLegendMaps(names, colors);
testCase.verifyError(@() CIMPQuantPlotter.plotXicGroupWithLayout( ...
    ric, total_xic, [10, 12], names, output_base, color_map, legend_map, layout), ...
    'CIMPQuantPlotter:LegendLayoutClipped');
for ext = {'.png', '.pdf', '.svg'}
    testCase.verifyFalse(isfile([output_base, ext{1}]), ...
        'Failed layout should not leave export artifacts.');
end
end

function testPreExportLegendLayoutFailureRemovesStaleArtifacts(testCase)
layout = CIMPQuantPlotter.getXicLegendLayoutConfig();
layout.figure_width_px = 180;
layout.figure_height_px = 180;
layout.left_margin_px = 80;
layout.bottom_margin_px = 50;
layout.top_margin_px = 40;
layout.axes_width_px = 100;
layout.axes_legend_gap_px = 40;
layout.legend_max_width_px = 20;
layout.legend_min_width_px = 20;
layout.legend_padding_px = 16;
layout.legend_min_line_chars = 4;

output_base = fullfile(testCase.TestData.output_dir, 'pre_export_layout_failure');
writeStaleXicArtifacts(testCase, output_base);

rt = linspace(10, 12, 40);
total_xic = {rt, exp(-((rt - 11).^2) / 0.08)};
ric = {rt, total_xic{2}};
names = {'WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW'};
colors = [0.0000, 0.4470, 0.7410];

[color_map, legend_map] = buildXicLegendMaps(names, colors);
testCase.verifyError(@() CIMPQuantPlotter.plotXicGroupWithLayout( ...
    ric, total_xic, [10, 12], names, output_base, color_map, legend_map, layout), ...
    'CIMPQuantPlotter:LegendLayoutClipped');
for ext = {'.png', '.pdf', '.svg'}
    testCase.verifyFalse(isfile([output_base, ext{1}]), ...
        'Pre-export layout failure should remove stale export artifacts.');
end
end

function testProductionLayoutExportStructure(testCase)
test_dir = fileparts(mfilename('fullpath'));
fixture_path = fullfile(test_dir, 'fixtures', 'xic_legend_layout', 'production_synthetic.mat');
testCase.assertTrue(isfile(fixture_path), ...
    'Missing fixture. Run fixtures/xic_legend_layout/create_production_synthetic_fixture.m first.');

fixture = load(fixture_path);
layout = CIMPQuantPlotter.getXicLegendLayoutConfig();
testCase.verifyEqual(layout.figure_width_px, 2000);
testCase.verifyEqual(layout.axes_width_fraction, 0.7);
testCase.verifyEqual(layout.axes_width_px, 1217);
testCase.verifyEqual(layout.figure_height_px, 800);
testCase.verifyEqual(layout.legend_max_width_px, 521);
testCase.verifyEqual(layout.all_font_size, 7);

output_base = fullfile(testCase.TestData.output_dir, 'production_synthetic_export');
for ext = {'.png', '.pdf', '.svg'}
    out_file = [output_base, ext{1}];
    if exist(out_file, 'file')
        delete(out_file);
    end
end

[color_map, legend_map] = buildXicLegendMaps(fixture.names, fixture.colors);
diagnostics = CIMPQuantPlotter.plotXicGroupWithLayout( ...
    fixture.ric, fixture.total_xic, fixture.categorized_intervals, ...
    fixture.names, output_base, color_map, legend_map, layout);

testCase.verifyTrue(diagnostics.is_legend_inside_figure, diagnostics.summary);
testCase.verifyTrue(diagnostics.is_text_inside_figure, diagnostics.summary);
testCase.verifyTrue(diagnostics.is_text_clear_of_axes, diagnostics.summary);
testCase.verifyGreaterThan(diagnostics.min_text_right_margin_px, 0, diagnostics.summary);
testCase.verifyGreaterThanOrEqual(diagnostics.rendered_legend_right_margin_px, 8, diagnostics.summary);

plot_band_center_y = diagnostics.axes_center_y_px;
text_center_y = diagnostics.legend_center_y_px;
testCase.verifyLessThan(abs(text_center_y - plot_band_center_y), 2, diagnostics.summary);

metrics = assertXicLegendExports(output_base, layout, diagnostics);
testCase.verifyGreaterThan(metrics.legend_right_margin_px, 0);
testCase.verifyGreaterThan(metrics.ink_row_fraction, 0.20);
testCase.verifyLessThan(metrics.ink_row_fraction, 0.85);
testCase.verifyLessThanOrEqual(metrics.right_blank_fraction, 0.20);
testCase.verifyTrue(metrics.svg_checked);
end

function writeStaleXicArtifacts(testCase, output_base)
for ext = {'.png', '.pdf', '.svg'}
    fid = fopen([output_base, ext{1}], 'w');
    testCase.assertGreaterThan(fid, 0, ['Could not create stale ', ext{1}, ' fixture.']);
    cleanup_obj = onCleanup(@() fclose(fid));
    fwrite(fid, 'stale export from previous run');
    clear cleanup_obj;
end
end
