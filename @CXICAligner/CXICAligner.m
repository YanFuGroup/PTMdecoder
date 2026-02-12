classdef CXICAligner < handle
    % Orchestrate anchor selection, transform fitting, and peak alignment

    properties (Access=private)
        m_anchor_selector
        m_options
    end

    methods
        function obj = CXICAligner(anchor_selector, options)
            % Construct an aligner with an anchor selector and default options
            % Input:
            %   anchor_selector (CXICAlignAnchorSelector, optional)
            %       Selector used to build anchors from FDR results
            %   options (struct, optional)
            %       Default alignment options
            % Output:
            %   obj (CXICAligner)
            %       Aligner instance
            if nargin < 1 || isempty(anchor_selector)
                anchor_selector = CXICAlignAnchorSelector();
            end
            if nargin < 2
                options = struct();
            end
            obj.m_anchor_selector = anchor_selector;
            obj.m_options = options;
        end

        [pair_models, report] = fitPairModels(obj, fdr_filtered_result_path, ms12DatasetIO, align_pairs, options)

        raw_names = getRawNamesFromReport(obj, quant_report)

        key = pairKey(obj, ref_raw, target_raw)

        [rt_peak_a, rt_peak_b, stats] = alignImpPair(obj, ms12DatasetIO, raw_name_a, imp_info_a, ...
            raw_name_b, imp_info_b, model, options, minMSMSnum, alpha, resFilterThres)
    end
end
