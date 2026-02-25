classdef CPTMdecoderWorkflowRunner < handle
    % Top-level workflow runner for PTMdecoder.

    properties (Access = private)
        m_config
        m_stage_executor_registry
    end

    methods
        function obj = CPTMdecoderWorkflowRunner(workflow_config)
            if nargin < 1 || isempty(workflow_config)
                error('CPTMdecoderWorkflowRunner:MissingConfig', ...
                    'workflow_config must be provided.');
            end
            obj.m_config = workflow_config;
            obj.m_stage_executor_registry = obj.buildStageExecutorRegistry();
        end

        function run(obj)
            for iStage = 1:numel(obj.m_config.stages)
                stage = obj.m_config.stages{iStage};
                if ~stage.enabled
                    continue;
                end
                obj.executeStage(stage);
            end
        end
    end

    methods (Access = private)
        function executeStage(obj, stage)
            if ~obj.m_stage_executor_registry.isKey(stage.name)
                error('CPTMdecoderWorkflowRunner:UnknownStage', ...
                    'Unknown workflow stage: %s', stage.name);
            end
            executor = obj.m_stage_executor_registry(stage.name);
            executor(stage);
        end

        function executors = buildStageExecutorRegistry(obj)
            executors = containers.Map('KeyType', 'char', 'ValueType', 'any');
            executors(CPTMdecoderWorkflowConfig.STAGE_MSMS_WORKFLOW) = @(stage) obj.executeMsmsStage(stage);
            executors(CPTMdecoderWorkflowConfig.STAGE_SITE_LEVEL) = @(stage) obj.executeSiteLevelStage(stage);
            executors(CPTMdecoderWorkflowConfig.STAGE_MERGE_TO_PAIR_LEVEL) = @(stage) obj.executeMergeToPairStage(stage);
            executors(CPTMdecoderWorkflowConfig.STAGE_MERGE_PAIRS_LEVEL) = @(stage) obj.executeMergePairsStage(stage);
        end

        function executeMsmsStage(~, stage)
            msms_cfg = stage.config;
            if isempty(msms_cfg)
                error('CPTMdecoderWorkflowRunner:MissingMsmsConfig', ...
                    'MSMS stage config is required.');
            end

            processor = CMSMSPepDeconv(msms_cfg);
            switch stage.action
                case CPTMdecoderWorkflowConfig.ACTION_MSMS_PEPTIDE
                    processor.runMsmsLevel();
                    processor.runImpQuantLevel();
                case CPTMdecoderWorkflowConfig.ACTION_MSMS_ONLY
                    processor.runMsmsLevel();
                case CPTMdecoderWorkflowConfig.ACTION_PEPTIDE_REQUANT
                    processor.runImpRequantLevel();
                case CPTMdecoderWorkflowConfig.ACTION_PEPTIDE_ONLY
                    processor.runImpQuantLevel();
                otherwise
                    error('CPTMdecoderWorkflowRunner:InvalidMsmsAction', ...
                        'Invalid msms workflow action: %s', stage.action);
            end
        end

        function executeSiteLevelStage(~, stage)
            site_cfg = stage.config;
            if isempty(site_cfg)
                error('CPTMdecoderWorkflowRunner:MissingSiteConfig', ...
                    'Site-level stage config is required.');
            end
            process = CSiteLevelSummary(site_cfg);
            process.run();
        end

        function executeMergeToPairStage(~, stage)
            pair_cfgs = stage.config;
            if isempty(pair_cfgs)
                return;
            end

            for idx_pairs = 1:numel(pair_cfgs)
                process = CMergeEachPair(pair_cfgs{idx_pairs});
                process.run();
            end
        end

        function executeMergePairsStage(~, stage)
            cfg = stage.config;
            if isempty(cfg)
                error('CPTMdecoderWorkflowRunner:MissingMergePairsConfig', ...
                    'Merge-pairs stage config is required.');
            end
            process = CMergePairs(cfg);
            process.run();
        end
    end
end
