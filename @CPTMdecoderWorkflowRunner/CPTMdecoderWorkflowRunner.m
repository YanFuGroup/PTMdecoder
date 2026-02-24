classdef CPTMdecoderWorkflowRunner < handle
    % Top-level workflow runner for PTMdecoder.

    properties (Access = private)
        m_config
        m_stage_executors
    end

    methods
        function obj = CPTMdecoderWorkflowRunner(workflow_config)
            if nargin < 1 || isempty(workflow_config)
                error('CPTMdecoderWorkflowRunner:MissingConfig', ...
                    'workflow_config must be provided.');
            end
            obj.m_config = workflow_config;
            obj.m_stage_executors = obj.buildStageExecutorMap();
        end

        function run(obj)
            for iStage = 1:numel(obj.m_config.stages)
                stage = obj.m_config.stages{iStage};
                if ~stage.enabled
                    continue;
                end
                obj.runStage(stage);
            end
        end
    end

    methods (Access = private)
        function runStage(obj, stage)
            if ~obj.m_stage_executors.isKey(stage.name)
                error('CPTMdecoderWorkflowRunner:UnknownStage', ...
                    'Unknown workflow stage: %s', stage.name);
            end
            executor = obj.m_stage_executors(stage.name);
            executor(stage);
        end

        function executors = buildStageExecutorMap(obj)
            executors = containers.Map('KeyType', 'char', 'ValueType', 'any');
            executors(CPTMdecoderWorkflowConfig.STAGE_MSMS_WORKFLOW) = @(stage) obj.runMsmsWorkflow(stage);
            executors(CPTMdecoderWorkflowConfig.STAGE_SITE_LEVEL) = @(stage) obj.runSiteLevelWorkflow(stage);
            executors(CPTMdecoderWorkflowConfig.STAGE_MERGE_TO_PAIR_LEVEL) = @(stage) obj.runMergeToPairWorkflow(stage);
            executors(CPTMdecoderWorkflowConfig.STAGE_MERGE_PAIRS_LEVEL) = @(stage) obj.runMergePairsWorkflow(stage);
        end

        function runMsmsWorkflow(~, stage)
            msms_cfg = stage.config;
            if isempty(msms_cfg)
                error('CPTMdecoderWorkflowRunner:MissingMsmsConfig', ...
                    'MSMS stage config is required.');
            end

            processor = CMSMSPepDeconv(msms_cfg);
            switch stage.action
                case CPTMdecoderWorkflowConfig.ACTION_MSMS_PEPTIDE
                    processor.startRun();
                case CPTMdecoderWorkflowConfig.ACTION_PEPTIDE_REQUANT
                    processor.requant_IMP();
                case CPTMdecoderWorkflowConfig.ACTION_PEPTIDE_ONLY
                    processor.quant_IMP();
                otherwise
                    error('CPTMdecoderWorkflowRunner:InvalidMsmsAction', ...
                        'Invalid msms workflow action: %s', stage.action);
            end
        end

        function runSiteLevelWorkflow(~, stage)
            site_cfg = stage.config;
            if isempty(site_cfg)
                error('CPTMdecoderWorkflowRunner:MissingSiteConfig', ...
                    'Site-level stage config is required.');
            end
            process = CSiteLevelSummary(site_cfg);
            process.summary_and_write();
        end

        function runMergeToPairWorkflow(~, stage)
            pair_cfgs = stage.config;
            if isempty(pair_cfgs)
                return;
            end

            for idx_pairs = 1:numel(pair_cfgs)
                process = CMergeEachPair(pair_cfgs{idx_pairs});
                process.merge_and_write();
            end
        end

        function runMergePairsWorkflow(~, stage)
            cfg = stage.config;
            if isempty(cfg)
                error('CPTMdecoderWorkflowRunner:MissingMergePairsConfig', ...
                    'Merge-pairs stage config is required.');
            end
            process = CMergePairs(cfg);
            process.merge_and_write();
        end
    end
end
