function [pep_rtrange_map, report] = buildAlignedRtrangeMap(obj, ...
    msms_result, fdr_filtered_result_path, buildRawIdentManagerFn)
% Build aligned RT range map using anchors and alignment models.
% Input:
%   obj (CIMPRequantAlignmentPipeline)
%       Pipeline instance
%   msms_result (CMSMSResult)
%       MSMS results from report_msms.txt
%   fdr_filtered_result_path (1 x 1 char/string)
%       FDR filtered result file path
%   buildRawIdentManagerFn (function_handle)
%       Function to build CIMPRawIdentManager from a spectrum_list
% Output:
%   pep_rtrange_map (containers.Map)
%       Map of record key -> aligned rt_peaks
%   report (struct)
%       Alignment stats and pair model info

if nargin < 4 || isempty(buildRawIdentManagerFn)
    error('buildAlignedRtrangeMap requires a raw identification builder function.');
end

raw_names = obj.getRawNamesFromMsmsResult(msms_result);
align_pairs = obj.m_align_strategy.getAlignmentPairs(raw_names);

[obj.m_pair_models, report] = obj.m_aligner.fitPairModels( ...
    fdr_filtered_result_path, obj.m_ms12DatasetIO, align_pairs, obj.m_align_options);

pep_rtrange_map = containers.Map();
state = struct('pep_rtrange_map', pep_rtrange_map, ...
    'num_groups', 0, 'num_aligned', 0, 'num_missing_raw_pair', 0);

if isempty(align_pairs)
    report.summary = state;
    return;
end

fprintf('Aligning RT ranges for %d run pairs...', size(align_pairs, 1));
print_progress = CPrintProgress(length(msms_result.Peptides));
for idx_psf = 1:length(msms_result.Peptides)
    print_progress = print_progress.update_show(idx_psf);
    rawIdentManager = buildRawIdentManagerFn(msms_result.Peptides(idx_psf).spectrum_list);
    raw_imps_map = obj.buildRawImpsMap(rawIdentManager);

    for idx_pair = 1:size(align_pairs, 1)
        raw_a = align_pairs{idx_pair, 1};
        raw_b = align_pairs{idx_pair, 2};
        if ~isKey(raw_imps_map, raw_a) || ~isKey(raw_imps_map, raw_b)
            state.num_missing_raw_pair = state.num_missing_raw_pair + 1;
            continue;
        end
        a_items = raw_imps_map(raw_a);
        b_items = raw_imps_map(raw_b);
        a_keys = a_items.keys;
        b_keys = b_items.keys;
        shared_keys = intersect(a_keys, b_keys);
        if isempty(shared_keys)
            state.num_missing_raw_pair = state.num_missing_raw_pair + 1;
            continue;
        end
        % shared_keys are IMP+charge observed in both raw files

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

            [rt_a, rt_b, stats] = obj.m_aligner.alignImpPair(...
                obj.m_ms12DatasetIO, raw_a, item_a, ...
                raw_b, item_b, model, ...
                obj.m_align_options, obj.m_minMSMSnum, obj.m_alpha, obj.m_resFilterThres);

            state.num_groups = state.num_groups + 1;
            state.num_aligned = state.num_aligned + stats.num_aligned;

            if ~isempty(rt_a)
                key = CIMPQuantRecord.build_id(item_a.impName, item_a.charge, raw_a);
                state.pep_rtrange_map(key) = rt_a;
            end
            if ~isempty(rt_b)
                key = CIMPQuantRecord.build_id(item_b.impName, item_b.charge, raw_b);
                state.pep_rtrange_map(key) = rt_b;
            end
        end
    end
end
print_progress.last_update();

pep_rtrange_map = state.pep_rtrange_map;
report.summary = state;
end
