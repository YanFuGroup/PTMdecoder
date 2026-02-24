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

            if isempty(obj.m_config.msms_peptide_config)
                error('CPTMdecoderWorkflowRunner:MissingMsmsConfig', ...
                    'msms_peptide_config is required when msms_workflow_mode is not none.');
            end

            processor = CMSMSPepDeconv(obj.m_config.msms_peptide_config);
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

            site_cfg = obj.m_config.site_level_config;
            if isempty(site_cfg)
                error('CPTMdecoderWorkflowRunner:MissingSiteConfig', ...
                    'site_level_config is required when site_level_on is true.');
            end
            process = CSiteLevelSummary(site_cfg);
            process.summary_and_write();
        end

        function runMergeToPairWorkflow(obj)
            if ~obj.m_config.merge_to_pair_level_on
                return;
            end

            pair_cfgs = obj.m_config.merge_to_pair_config;
            if isempty(pair_cfgs)
                return;
            end

            for idx_pairs = 1:numel(pair_cfgs)
                process = CMergeEachPair(pair_cfgs{idx_pairs});
                process.merge_and_write();
            end
        end

        function runMergePairsWorkflow(obj)
            if ~obj.m_config.merge_pairs_level_on
                return;
            end

            cfg = obj.m_config.merge_pairs_config;
            if isempty(cfg)
                error('CPTMdecoderWorkflowRunner:MissingMergePairsConfig', ...
                    'merge_pairs_config is required when merge_pairs_level_on is true.');
            end
            process = CMergePairs(cfg);
            process.merge_and_write();
        end
    end
end
