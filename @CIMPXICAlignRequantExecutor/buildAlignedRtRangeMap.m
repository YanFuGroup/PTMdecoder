function [pep_rtrange_map, report] = buildAlignedRtRangeMap(obj, ...
    fdr_filtered_result_path, rawIdentManagers, base_pep_rtrange_map, base_groups_by_peptide)
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
%   base_groups_by_peptide (1 x N cell, optional)
%       Prebuilt CIMPGroup arrays for each peptide.
%       When provided, align executor reuses grouped contexts and skips
%       building raw IMP maps from rawIdentManagers again.
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
if nargin < 5
    base_groups_by_peptide = [];
end
if ~isa(base_pep_rtrange_map, 'containers.Map')
    error('buildAlignedRtRangeMap requires base_pep_rtrange_map to be a containers.Map.');
end
if ~isempty(base_groups_by_peptide) && numel(base_groups_by_peptide) ~= numel(rawIdentManagers)
    error('CIMPXICAlignRequantExecutor:InvalidBaseGroupsByPeptide', ...
        'base_groups_by_peptide must match rawIdentManagers in length.');
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

CLogger.info('Aligning XIC candidates for %d run pairs...', size(align_pairs, 1));
print_progress = CPrintProgress(length(rawIdentManagers), 'align_xic_candidates');
for idx_psf = 1:length(rawIdentManagers)
    print_progress = print_progress.update_show(idx_psf);
    if ~isempty(base_groups_by_peptide)
        raw_imps_map = obj.buildRawImpsMapFromBaseGroups(base_groups_by_peptide{idx_psf});
    else
        rawIdentManager = rawIdentManagers{idx_psf};
        raw_imps_map = obj.buildRawImpsMap(rawIdentManager);
    end

    for idx_pair = 1:size(align_pairs, 1)
        raw_a = align_pairs{idx_pair, 1};
        raw_b = align_pairs{idx_pair, 2};

        if ~isKey(raw_imps_map, raw_a) || ~isKey(raw_imps_map, raw_b)
            if isKey(raw_imps_map, raw_a)
                a_items = raw_imps_map(raw_a);
                [state.pep_rtrange_map, applied_a] = applyUnalignedForBaseKeys(...
                    state.pep_rtrange_map, collectBaseKeys(a_items), unaligned_action);
                state.num_unaligned_action_applied = state.num_unaligned_action_applied + applied_a;
            end
            if isKey(raw_imps_map, raw_b)
                b_items = raw_imps_map(raw_b);
                [state.pep_rtrange_map, applied_b] = applyUnalignedForBaseKeys(...
                    state.pep_rtrange_map, collectBaseKeys(b_items), unaligned_action);
                state.num_unaligned_action_applied = state.num_unaligned_action_applied + applied_b;
            end
            state.num_missing_raw_pair = state.num_missing_raw_pair + 1;
            continue;
        end
        a_items = raw_imps_map(raw_a);
        b_items = raw_imps_map(raw_b);
        [shared_entries, unaligned_a_keys, unaligned_b_keys] = ...
            buildPairProcessingPlan(a_items, b_items, base_pep_rtrange_map);

        [state.pep_rtrange_map, applied_a] = applyUnalignedForBaseKeys(...
            state.pep_rtrange_map, unaligned_a_keys, unaligned_action);
        state.num_unaligned_action_applied = state.num_unaligned_action_applied + applied_a;
        [state.pep_rtrange_map, applied_b] = applyUnalignedForBaseKeys(...
            state.pep_rtrange_map, unaligned_b_keys, unaligned_action);
        state.num_unaligned_action_applied = state.num_unaligned_action_applied + applied_b;

        if isempty(shared_entries)
            state.num_missing_raw_pair = state.num_missing_raw_pair + 1;
            continue;
        end

        model_key = obj.m_aligner.pairKey(raw_a, raw_b);
        if ~obj.m_pair_models.isKey(model_key)
            CLogger.warn('No alignment model for pair (%s, %s). Skipping group alignment.', raw_a, raw_b);
            continue;
        end
        model = obj.m_pair_models(model_key);

        for idx_key = 1:numel(shared_entries)
            shared_entry = shared_entries(idx_key);
            item_a = shared_entry.item_a;
            item_b = shared_entry.item_b;
            key_a = shared_entry.key_a;
            key_b = shared_entry.key_b;
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
CLogger.info('Aligning XIC candidates done.');

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

function [pep_rtrange_map, applied_total] = applyUnalignedForBaseKeys(pep_rtrange_map, base_keys, unaligned_action)
applied_total = 0;
for idx_key = 1:numel(base_keys)
    key = base_keys{idx_key};
    [pep_rtrange_map, applied] = applyUnalignedAction(pep_rtrange_map, key, unaligned_action);
    applied_total = applied_total + applied;
end
end


function base_keys = collectBaseKeys(item_map)
item_keys = item_map.keys;
base_keys = cell(numel(item_keys), 1);
for idx_key = 1:numel(item_keys)
    item = item_map(item_keys{idx_key});
    base_keys{idx_key} = getItemBaseKey(item);
end
end

function [shared_entries, unaligned_a_keys, unaligned_b_keys] = ...
    buildPairProcessingPlan(a_items, b_items, base_pep_rtrange_map)
shared_entries = repmat(struct('item_a', [], 'item_b', [], 'key_a', '', 'key_b', ''), 0, 1);
unaligned_a_keys = cell(0, 1);
unaligned_b_keys = cell(0, 1);

a_keys = a_items.keys;
b_keys = b_items.keys;

for idx_key = 1:numel(a_keys)
    imp_key = a_keys{idx_key};
    item_a = a_items(imp_key);
    key_a = getItemBaseKey(item_a);

    if ~isKey(b_items, imp_key)
        unaligned_a_keys{end + 1, 1} = key_a; %#ok<AGROW>
        continue;
    end

    item_b = b_items(imp_key);
    key_b = getItemBaseKey(item_b);
    if base_pep_rtrange_map.isKey(key_a) && base_pep_rtrange_map.isKey(key_b)
        shared_entries(end + 1, 1) = struct( ...
            'item_a', item_a, ...
            'item_b', item_b, ...
            'key_a', key_a, ...
            'key_b', key_b); %#ok<AGROW>
    else
        unaligned_a_keys{end + 1, 1} = key_a; %#ok<AGROW>
        unaligned_b_keys{end + 1, 1} = key_b; %#ok<AGROW>
    end
end

for idx_key = 1:numel(b_keys)
    imp_key = b_keys{idx_key};
    if isKey(a_items, imp_key)
        continue;
    end
    item_b = b_items(imp_key);
    unaligned_b_keys{end + 1, 1} = getItemBaseKey(item_b); %#ok<AGROW>
end
end


function base_key = getItemBaseKey(item)
if ~isfield(item, 'baseKey') || isempty(item.baseKey)
    error('CIMPXICAlignRequantExecutor:MissingItemBaseKey', ...
        'imp_info must contain precomputed baseKey for alignment processing.');
end
base_key = item.baseKey;
end
