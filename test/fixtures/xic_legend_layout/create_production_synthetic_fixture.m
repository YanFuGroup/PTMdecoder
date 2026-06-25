% One-time script to create production_synthetic.mat fixture for legend layout tests.
fixture_dir = fileparts(mfilename('fullpath'));

rt = linspace(10, 12, 120);
total_xic = {rt, exp(-((rt - 11).^2) / 0.06)};
ric = {
    rt, 0.55 .* total_xic{2};
    rt, 0.35 .* total_xic{2}
};
categorized_intervals = [10, 12];
names = {
    '_PEPT{phospho}IDESK_WITH_EXTRA_LONG_CONTEXT_FOR_LAYOUT_';
    '_PEPTIDES{phospho}K_WITH_EXTRA_LONG_CONTEXT_FOR_LAYOUT_'
};
colors = [
    0.0000, 0.4470, 0.7410;
    0.8500, 0.3250, 0.0980
];

fixture_path = fullfile(fixture_dir, 'production_synthetic.mat');
save(fixture_path, 'ric', 'total_xic', 'categorized_intervals', 'names', 'colors');
fprintf('Wrote fixture: %s\n', fixture_path);
