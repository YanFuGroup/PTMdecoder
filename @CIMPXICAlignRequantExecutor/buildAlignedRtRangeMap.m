function [pep_rtrange_map, report] = buildAlignedRtRangeMap(obj, ...
    fdr_filtered_result_path, rawIdentManagers, base_pep_rtrange_map)
% Build aligned RT range map using alignment models on top of base quant ranges.
%
% Function role:
%   1) Fit run-pair alignment models from FDR-filtered anchors.
%   2) For each shared IMP across aligned run pairs, evaluate candidate RT peaks
%      constrained by base_pep_rtrange_map values.
%   3) Update pep_rtrange_map with aligned peak selections.
%   4) For IMPs without valid aligned peak pairs, apply configurable policy:
%        - 'drop'      : remove IMP key from pep_rtrange_map (default)
%        - 'zero-label': keep IMP key but set all check_label = 0
%
% Input:
%   obj (CIMPXICAlignRequantExecutor)
%       Executor instance
%   fdr_filtered_result_path (1 x 1 char/string)
%       FDR filtered result file path
%   rawIdentManagers (1 x N cell)
%       Per-peptide CIMPRawIdentManager list (used to build per-raw IMP info)
%   base_pep_rtrange_map (containers.Map)
%       Base peptide RT map from initial quantification.
%       Key format: CIMPQuantRecord.build_id(imp_name, charge, raw_name)
%       Value: struct array with fields rt_start/rt_end/ratio/check_label.
%       This map defines the candidate RT interval universe for alignment.
% Output:
%   pep_rtrange_map (containers.Map)
%       Updated RT map after alignment and unaligned policy handling.
%   report (struct)
%       Alignment report including fitted pair models and summary counters:
%         - num_groups
%         - num_aligned
%         - num_missing_raw_pair
%         - num_updated_from_base
%         - num_unaligned_action_applied
%         - unaligned_imp_action

if nargin < 3 || isempty(rawIdentManagers)
    error('buildAlignedRtRangeMap requires prebuilt rawIdentManagers.');
end
if nargin < 4 || isempty(base_pep_rtrange_map)
    error('buildAlignedRtRangeMap requires base_pep_rtrange_map from initial quantification.');
end
if ~isa(base_pep_rtrange_map, 'containers.Map')
    error('buildAlignedRtRangeMap requires base_pep_rtrange_map to be a containers.Map.');
end

raw_names = obj.getRawNamesFromRawIdentManagers(rawIdentManagers);
align_pairs = obj.m_align_strategy.getAlignmentPairs(raw_names);

[obj.m_pair_models, report] = obj.m_aligner.fitPairModels( ...
    fdr_filtered_result_path, obj.m_ms12DatasetIO, align_pairs, obj.m_align_options);

unaligned_action = lower(strtrim(CStructOptionUtils.get(obj.m_align_options, ...
    'unaligned_imp_action', 'drop')));
if ~any(strcmp(unaligned_action, {'drop', 'zero-label'}))
    error('CIMPXICAlignRequantExecutor:InvalidUnalignedImpAction', ...
        'align_options.unaligned_imp_action must be ''drop'' or ''zero-label'', got: %s', unaligned_action);
end

pep_rtrange_map = copyRtRangeMap(base_pep_rtrange_map);
state = struct('pep_rtrange_map', pep_rtrange_map, ...
    'num_groups', 0, 'num_aligned', 0, 'num_missing_raw_pair', 0, ...
    'num_updated_from_base', 0, 'num_unaligned_action_applied', 0, ...
    'unaligned_imp_action', unaligned_action);

if isempty(align_pairs)
    report.summary = state;
    return;
end

fprintf('Aligning XIC candidates for %d run pairs...', size(align_pairs, 1));
print_progress = CPrintProgress(length(rawIdentManagers));
for idx_psf = 1:length(rawIdentManagers)
    print_progress = print_progress.update_show(idx_psf);
    rawIdentManager = rawIdentManagers{idx_psf};
    raw_imps_map = obj.buildRawImpsMap(rawIdentManager);

    for idx_pair = 1:size(align_pairs, 1)
        raw_a = align_pairs{idx_pair, 1};
        raw_b = align_pairs{idx_pair, 2};

        if ~isKey(raw_imps_map, raw_a) || ~isKey(raw_imps_map, raw_b)
            if isKey(raw_imps_map, raw_a)
                a_items = raw_imps_map(raw_a);
                [state.pep_rtrange_map, applied_a] = applyUnalignedForImpKeys(...
                    state.pep_rtrange_map, a_items.keys, raw_a, unaligned_action);
                state.num_unaligned_action_applied = state.num_unaligned_action_applied + applied_a;
            end
            if isKey(raw_imps_map, raw_b)
                b_items = raw_imps_map(raw_b);
                [state.pep_rtrange_map, applied_b] = applyUnalignedForImpKeys(...
                    state.pep_rtrange_map, b_items.keys, raw_b, unaligned_action);
                state.num_unaligned_action_applied = state.num_unaligned_action_applied + applied_b;
            end
            state.num_missing_raw_pair = state.num_missing_raw_pair + 1;
            continue;
        end
        a_items = raw_imps_map(raw_a);
        b_items = raw_imps_map(raw_b);
        a_keys = a_items.keys;
        b_keys = b_items.keys;
        shared_keys = intersect(a_keys, b_keys);
        shared_keys = filterSharedKeysByBaseMap(shared_keys, raw_a, raw_b, base_pep_rtrange_map);

        [state.pep_rtrange_map, applied_a] = applyUnalignedForImpKeys(...
            state.pep_rtrange_map, setdiff(a_keys, shared_keys), raw_a, unaligned_action);
        state.num_unaligned_action_applied = state.num_unaligned_action_applied + applied_a;
        [state.pep_rtrange_map, applied_b] = applyUnalignedForImpKeys(...
            state.pep_rtrange_map, setdiff(b_keys, shared_keys), raw_b, unaligned_action);
        state.num_unaligned_action_applied = state.num_unaligned_action_applied + applied_b;

        if isempty(shared_keys)
            state.num_missing_raw_pair = state.num_missing_raw_pair + 1;
            continue;
        end

        model_key = obj.m_aligner.pairKey(raw_a, raw_b);
        if ~obj.m_pair_models.isKey(model_key)
            warning('No alignment model for pair (%s, %s). Skipping group alignment.', raw_a, raw_b);
            continue;
        end
        model = obj.m_pair_models(model_key);

        for idx_key = 1:numel(shared_keys)
            imp_key = shared_keys{idx_key};
            item_a = a_items(imp_key);
            item_b = b_items(imp_key);
            key_a = CIMPQuantRecord.build_id(item_a.impName, item_a.charge, raw_a);
            key_b = CIMPQuantRecord.build_id(item_b.impName, item_b.charge, raw_b);
            item_a.candidateRtPeaks = base_pep_rtrange_map(key_a);
            item_b.candidateRtPeaks = base_pep_rtrange_map(key_b);

            [rt_a, rt_b, stats] = obj.m_aligner.alignImpPair(...
                obj.m_ms12DatasetIO, raw_a, item_a, ...
                raw_b, item_b, model, ...
                obj.m_align_options, obj.m_minMSMSnum, obj.m_resFilterThres);

            state.num_groups = state.num_groups + 1;
            state.num_aligned = state.num_aligned + stats.num_aligned;

            if ~isempty(rt_a)
                if state.pep_rtrange_map.isKey(key_a)
                    state.pep_rtrange_map(key_a) = rt_a;
                    state.num_updated_from_base = state.num_updated_from_base + 1;
                end
            else
                [state.pep_rtrange_map, applied] = applyUnalignedAction(...
                    state.pep_rtrange_map, key_a, unaligned_action);
                state.num_unaligned_action_applied = state.num_unaligned_action_applied + applied;
            end
            if ~isempty(rt_b)
                if state.pep_rtrange_map.isKey(key_b)
                    state.pep_rtrange_map(key_b) = rt_b;
                    state.num_updated_from_base = state.num_updated_from_base + 1;
                end
            else
                [state.pep_rtrange_map, applied] = applyUnalignedAction(...
                    state.pep_rtrange_map, key_b, unaligned_action);
                state.num_unaligned_action_applied = state.num_unaligned_action_applied + applied;
            end
        end
    end
end
print_progress.last_update();
fprintf('done.\n');

pep_rtrange_map = state.pep_rtrange_map;
report.summary = state;
end

function cloned_map = copyRtRangeMap(source_map)
cloned_map = containers.Map();
source_keys = source_map.keys;
for idx = 1:numel(source_keys)
    key = source_keys{idx};
    cloned_map(key) = source_map(key);
end
end

function [pep_rtrange_map, applied] = applyUnalignedAction(pep_rtrange_map, key, unaligned_action)
% Apply policy for IMPs without valid aligned peak selection.
% Input:
%   pep_rtrange_map (containers.Map)
%       Current output map
%   key (1 x 1 char/string)
%       IMP key in CIMPQuantRecord.build_id format
%   unaligned_action (1 x 1 char/string)
%       'drop' or 'zero-label'
% Output:
%   pep_rtrange_map (containers.Map)
%       Updated output map
%   applied (1 x 1 double/int)
%       1 if action applied, 0 otherwise
applied = 0;
if ~pep_rtrange_map.isKey(key)
    return;
end

switch unaligned_action
    case 'drop'
        remove(pep_rtrange_map, key);
        applied = 1;
    case 'zero-label'
        peaks = pep_rtrange_map(key);
        if isempty(peaks)
            return;
        end
        for idx_peak = 1:numel(peaks)
            peaks(idx_peak).check_label = 0;
        end
        pep_rtrange_map(key) = peaks;
        applied = 1;
end
end

function [pep_rtrange_map, applied_total] = applyUnalignedForImpKeys(pep_rtrange_map, imp_keys, raw_name, unaligned_action)
applied_total = 0;
for idx_key = 1:numel(imp_keys)
    [imp_name, imp_charge, ok] = parseImpKey(imp_keys{idx_key});
    if ~ok
        continue;
    end
    key = CIMPQuantRecord.build_id(imp_name, imp_charge, raw_name);
    [pep_rtrange_map, applied] = applyUnalignedAction(pep_rtrange_map, key, unaligned_action);
    applied_total = applied_total + applied;
end
end

function [imp_name, imp_charge, is_ok] = parseImpKey(imp_key)
imp_name = '';
imp_charge = NaN;
is_ok = false;

split_keys = strsplit(imp_key, '|');
if numel(split_keys) < 2
    return;
end
imp_name = split_keys{1};
imp_charge = str2double(split_keys{2});
if isnan(imp_charge)
    return;
end
is_ok = true;
end

function filtered_shared_keys = filterSharedKeysByBaseMap(shared_keys, raw_a, raw_b, base_pep_rtrange_map)
% Filter shared IMP keys to those present in base map for both runs.
% Input:
%   shared_keys (cell)
%       Shared imp keys in 'impName|charge' format
%   raw_a/raw_b (1 x 1 char/string)
%       Run names of current pair
%   base_pep_rtrange_map (containers.Map)
%       Base RT range map
% Output:
%   filtered_shared_keys (cell)
%       Keys with both raw_a and raw_b entries in base map
filtered_shared_keys = cell(0, 1);
for idx_key = 1:numel(shared_keys)
    imp_key = shared_keys{idx_key};
    split_keys = strsplit(imp_key, '|');
    if numel(split_keys) < 2
        continue;
    end
    imp_name = split_keys{1};
    imp_charge = str2double(split_keys{2});
    if isnan(imp_charge)
        continue;
    end

    key_a = CIMPQuantRecord.build_id(imp_name, imp_charge, raw_a);
    key_b = CIMPQuantRecord.build_id(imp_name, imp_charge, raw_b);
    if base_pep_rtrange_map.isKey(key_a) && base_pep_rtrange_map.isKey(key_b)
        filtered_shared_keys{end+1, 1} = imp_key; %#ok<AGROW>
    end
end
end
