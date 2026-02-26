classdef CPTMdecoderWorkflowRunner < handle
    % Top-level workflow runner for PTMdecoder.

    properties (Access = private)
        m_config
        m_stage_executor_registry
    end

    methods
        function obj = CPTMdecoderWorkflowRunner(workflow_config)
            % Input:
            %   workflow_config (CPTMdecoderWorkflowConfig)
            %       parsed workflow config object
            if nargin < 1 || isempty(workflow_config)
                error('CPTMdecoderWorkflowRunner:MissingConfig', ...
                    'workflow_config must be provided.');
            end
            obj.m_config = workflow_config;
            obj.m_stage_executor_registry = obj.buildStageExecutorRegistry();
        end

        function run(obj)
            % Run all enabled workflow stages in order.
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
            % Execute one workflow stage using registered executor.
            % Input:
            %   stage (struct)
            %       stage struct with fields: name, action, config, enabled
            %       - stage.name is used as runtime dispatch key
            %       - stage.enabled is checked by caller run()
            %       - stage.action is retained for compatibility tracing
            %       - stage.config is consumed by concrete stage executor
            if ~obj.m_stage_executor_registry.isKey(stage.name)
                error('CPTMdecoderWorkflowRunner:UnknownStage', ...
                    'Unknown workflow stage: %s', stage.name);
            end
            executor = obj.m_stage_executor_registry(stage.name);
            executor(stage);
        end

        function executors = buildStageExecutorRegistry(obj)
            % Build stage-name to executor-function registry.
            % Output:
            %   executors (containers.Map)
            %       map from stage name to execution function handle
            executors = containers.Map('KeyType', 'char', 'ValueType', 'any');
            executors(CPTMdecoderWorkflowConfig.STAGE_MSMS_LEVEL) = @(stage) obj.executeMsmsLevelStage(stage);
            executors(CPTMdecoderWorkflowConfig.STAGE_PEPTIDE_QUANT) = @(stage) obj.executePeptideQuantStage(stage);
            executors(CPTMdecoderWorkflowConfig.STAGE_PEPTIDE_REQUANT) = @(stage) obj.executePeptideRequantStage(stage);
            executors(CPTMdecoderWorkflowConfig.STAGE_NORM_PEPTIDE_QUANT) = @(stage) obj.executeNormPeptideQuantStage(stage);
            executors(CPTMdecoderWorkflowConfig.STAGE_NORM_PEPTIDE_REQUANT) = @(stage) obj.executeNormPeptideRequantStage(stage);
            executors(CPTMdecoderWorkflowConfig.STAGE_SITE_LEVEL) = @(stage) obj.executeSiteLevelStage(stage);
            executors(CPTMdecoderWorkflowConfig.STAGE_MERGE_TO_PAIR_LEVEL) = @(stage) obj.executeMergeToPairStage(stage);
            executors(CPTMdecoderWorkflowConfig.STAGE_MERGE_PAIRS_LEVEL) = @(stage) obj.executeMergePairsStage(stage);
        end

        function executeMsmsLevelStage(~, stage)
            % Execute MSMS-level stage.
            % Input:
            %   stage (struct)
            %       stage struct
            %       - stage.name should be STAGE_MSMS_LEVEL
            %       - stage.action should be ACTION_MSMS_ONLY
            %       - stage.config is MSMS config struct
            msms_cfg = stage.config;
            if isempty(msms_cfg)
                error('CPTMdecoderWorkflowRunner:MissingMsmsConfig', ...
                    'MSMS level stage config is required.');
            end

            service = CMSMSLevelService(msms_cfg);
            service.run();
        end

        function executePeptideQuantStage(~, stage)
            % Execute peptide quantification stage.
            % Input:
            %   stage (struct)
            %       stage struct
            %       - stage.name should be STAGE_PEPTIDE_QUANT
            %       - stage.action should be ACTION_PEPTIDE_ONLY
            %       - stage.config is MSMS config struct
            msms_cfg = stage.config;
            if isempty(msms_cfg)
                error('CPTMdecoderWorkflowRunner:MissingMsmsConfig', ...
                    'Peptide quant stage config is required.');
            end

            service = CPeptideQuantService(msms_cfg);
            service.run();
        end

        function executePeptideRequantStage(~, stage)
            % Execute peptide re-quantification stage.
            % Input:
            %   stage (struct)
            %       stage struct
            %       - stage.name should be STAGE_PEPTIDE_REQUANT
            %       - stage.action should be ACTION_PEPTIDE_REQUANT
            %       - stage.config is MSMS config struct
            msms_cfg = stage.config;
            if isempty(msms_cfg)
                error('CPTMdecoderWorkflowRunner:MissingMsmsConfig', ...
                    'Peptide requant stage config is required.');
            end

            service = CPeptideRequantService(msms_cfg);
            service.run();
        end

        function executeNormPeptideQuantStage(~, stage)
            % Execute normalization peptide quantification stage.
            % Input:
            %   stage (struct)
            %       stage struct
            %       - stage.name should be STAGE_NORM_PEPTIDE_QUANT
            %       - stage.action should be ACTION_NORM_PEPTIDE_QUANT
            %       - stage.config is normalization quant config struct
            cfg = stage.config;
            if isempty(cfg)
                error('CPTMdecoderWorkflowRunner:MissingNormQuantConfig', ...
                    'Normalization quant stage config is required.');
            end

            service = CNormalizationQuantService(cfg);
            service.run();
        end

        function executeNormPeptideRequantStage(~, stage)
            % Execute normalization peptide re-quantification stage.
            % Input:
            %   stage (struct)
            %       stage struct
            %       - stage.name should be STAGE_NORM_PEPTIDE_REQUANT
            %       - stage.action should be ACTION_NORM_PEPTIDE_REQUANT
            %       - stage.config is normalization requant config struct
            cfg = stage.config;
            if isempty(cfg)
                error('CPTMdecoderWorkflowRunner:MissingNormRequantConfig', ...
                    'Normalization requant stage config is required.');
            end

            service = CNormalizationRequantService(cfg);
            service.run();
        end

        function executeSiteLevelStage(~, stage)
            % Execute site-level summary stage.
            % Input:
            %   stage (struct)
            %       stage struct
            %       - stage.name should be STAGE_SITE_LEVEL
            %       - stage.action should be ACTION_SITE_SUMMARY
            %       - stage.config is site-level summary config struct
            site_cfg = stage.config;
            if isempty(site_cfg)
                error('CPTMdecoderWorkflowRunner:MissingSiteConfig', ...
                    'Site-level stage config is required.');
            end
            process = CSiteLevelSummary(site_cfg);
            process.run();
        end

        function executeMergeToPairStage(~, stage)
            % Execute merge-to-pair stage for each pair config.
            % Input:
            %   stage (struct)
            %       stage struct
            %       - stage.name should be STAGE_MERGE_TO_PAIR_LEVEL
            %       - stage.action should be ACTION_MERGE_EACH_PAIR
            %       - stage.config is pair config cell array
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
            % Execute merge-pairs stage.
            % Input:
            %   stage (struct)
            %       stage struct
            %       - stage.name should be STAGE_MERGE_PAIRS_LEVEL
            %       - stage.action should be ACTION_MERGE_PAIRS
            %       - stage.config is merge-pairs config struct
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
