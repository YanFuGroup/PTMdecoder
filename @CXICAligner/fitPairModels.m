function [pair_models, report] = fitPairModels(obj, fdr_filtered_result_path, ms12DatasetIO, align_pairs, options)
% Fit alignment models for all requested run pairs.
% Input:
%   obj (CXICAligner)
%       Aligner instance
%   fdr_filtered_result_path (1 x 1 char/string)
%       Path to FDR filtered result file
%   ms12DatasetIO (CMS12DatasetIO)
%       MS1/MS2 dataset IO (for RT resolution)
%   align_pairs (N x 2 cell)
%       [ref_raw, target_raw] pairs
%   options (struct, optional)
%       Alignment options (overrides defaults):
%         - min_psm (1 x 1 double)
%             Minimum PSM count per peptide anchor. Passed to CXICAlignAnchorSelector.
%             Default: 3.
%         - num_bins (1 x 1 double)
%             Number of bins for local offset fitting. Passed to CXICAlignTransform.fit.
%             Default: 5.
%         - min_per_bin (1 x 1 double)
%             Minimum anchors per bin for local offsets. Passed to CXICAlignTransform.fit.
%             Default: 5.
% Output:
%   pair_models (containers.Map)
%       Map of pairKey -> transform model
%   report (struct)
%       Pair stats and total anchor count

if nargin < 5
    options = obj.m_options;
end

anchors = obj.m_anchor_selector.selectAnchors(fdr_filtered_result_path, ms12DatasetIO, options);
pair_models = containers.Map('KeyType', 'char', 'ValueType', 'any');
report = struct('pairs', [], 'total_anchors', numel(anchors));

if isempty(align_pairs)
    return;
end

% Build peptide -> RT map per raw
raw_maps = containers.Map('KeyType', 'char', 'ValueType', 'any');
for idx = 1:numel(anchors)
    anchor = anchors(idx);
    if ~isKey(raw_maps, anchor.raw_name)
        raw_maps(anchor.raw_name) = containers.Map('KeyType', 'char', 'ValueType', 'double');
    end
    raw_map = raw_maps(anchor.raw_name);
    raw_map(anchor.peptide) = anchor.rt;
    raw_maps(anchor.raw_name) = raw_map;
end

pair_stats = repmat(struct('ref_raw', '', 'target_raw', '', 'num_anchors', 0, ...
    'slope', 1, 'intercept', 0, 'median_abs_resid', NaN), 0, 1);

for idx = 1:size(align_pairs, 1)
    ref_raw = align_pairs{idx, 1};
    target_raw = align_pairs{idx, 2};
    key = obj.pairKey(ref_raw, target_raw);

    if ~isKey(raw_maps, ref_raw) || ~isKey(raw_maps, target_raw)
        pair_models(key) = CXICAlignTransform.fit([], [], options);
        pair_stats(end+1) = struct('ref_raw', ref_raw, 'target_raw', target_raw, ...
            'num_anchors', 0, 'slope', 1, 'intercept', 0, 'median_abs_resid', NaN); %#ok<AGROW>
        continue;
    end

    % Get shared peptides and RTs for the pair
    ref_map = raw_maps(ref_raw);
    target_map = raw_maps(target_raw);
    ref_peps = ref_map.keys;
    target_peps = target_map.keys;
    shared_peps = intersect(ref_peps, target_peps);

    ref_rts = zeros(numel(shared_peps), 1);
    target_rts = zeros(numel(shared_peps), 1);
    for j = 1:numel(shared_peps)
        pep = shared_peps{j};
        ref_rts(j) = ref_map(pep);
        target_rts(j) = target_map(pep);
    end

    % Fit model for the pair
    model = CXICAlignTransform.fit(ref_rts, target_rts, options);
    pair_models(key) = model;

    % Calculate pair stats metrics for the report
    if model.has_model && ~isempty(shared_peps)
        pred = CXICAlignTransform.apply(model, ref_rts);
        resid = abs(target_rts - pred);
        med_resid = median(resid);
    else
        med_resid = NaN;
    end

    pair_stats(end+1) = struct('ref_raw', ref_raw, 'target_raw', target_raw, ...
        'num_anchors', numel(shared_peps), 'slope', model.slope, ...
        'intercept', model.intercept, 'median_abs_resid', med_resid);   %#ok<AGROW>
end

report.pairs = pair_stats;
end
