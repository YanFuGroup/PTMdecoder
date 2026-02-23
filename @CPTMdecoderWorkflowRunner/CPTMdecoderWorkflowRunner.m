classdef CPTMdecoderWorkflowRunner < handle
    % Top-level workflow runner for PTMdecoder.

    properties (Access = private)
        m_config
    end

    methods
        function obj = CPTMdecoderWorkflowRunner(workflow_config)
            if nargin < 1 || isempty(workflow_config)
                error('CPTMdecoderWorkflowRunner:MissingConfig', ...
                    'workflow_config must be provided.');
            end
            obj.m_config = workflow_config;
        end

        function run(obj)
            for iStage = 1:length(obj.m_config.execution_order)
                stage = obj.m_config.execution_order{iStage};
                switch stage
                    case 'msms_workflow'
                        obj.runMsmsWorkflow();
                    case 'site_level'
                        obj.runSiteLevelWorkflow();
                    case 'merge_to_pair_level'
                        obj.runMergeToPairWorkflow();
                    case 'merge_pairs_level'
                        obj.runMergePairsWorkflow();
                    otherwise
                        error('CPTMdecoderWorkflowRunner:UnknownStage', ...
                            'Unknown workflow stage: %s', stage);
                end
            end
        end
    end

    methods (Access = private)
        function runMsmsWorkflow(obj)
            mode = obj.m_config.msms_workflow_mode;
            if strcmp(mode, 'none')
                return;
            end

            processor = CMSMSPepDeconv(obj.m_config.legacy_task_param);
            switch mode
                case 'msms_peptide'
                    processor.startRun();
                case 'peptide_requant'
                    processor.requant_IMP();
                case 'peptide_only'
                    processor.quant_IMP();
                otherwise
                    error('CPTMdecoderWorkflowRunner:InvalidMsmsMode', ...
                        'Invalid msms workflow mode: %s', mode);
            end
        end

        function runSiteLevelWorkflow(obj)
            if ~obj.m_config.site_level_on
                return;
            end

            cfg = obj.m_config.site_level_config;
            process = CSiteLevelSummary( ...
                cfg.input_path, ...
                cfg.output_intere_path, ...
                cfg.output_unintere_path, ...
                cfg.protein_name_abbr, ...
                cfg.mod_name_abbr, ...
                cfg.ignore_strings, ...
                cfg.column_idxs);
            process.summary_and_write();
        end

        function runMergeToPairWorkflow(obj)
            if ~obj.m_config.merge_to_pair_level_on
                return;
            end

            cfg = obj.m_config.merge_to_pair_config;
            if isempty(cfg.pairs)
                return;
            end

            for idx_pairs = 1:size(cfg.pairs, 1)
                current_pair = cfg.pairs(idx_pairs, :);
                process = CMergeEachPair( ...
                    current_pair{1}, ...
                    current_pair{2}, ...
                    current_pair{3}, ...
                    {cfg.left_name, cfg.right_name}, ...
                    cfg.ignore_strings);
                process.merge_and_write();
            end
        end

        function runMergePairsWorkflow(obj)
            if ~obj.m_config.merge_pairs_level_on
                return;
            end

            cfg = obj.m_config.merge_pairs_config;
            process = CMergePairs(cfg.result_paths, cfg.output_path, cfg.group_titles);
            process.merge_and_write();
        end
    end
end
