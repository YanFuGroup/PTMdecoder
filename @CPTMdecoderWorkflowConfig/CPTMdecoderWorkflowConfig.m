classdef CPTMdecoderWorkflowConfig
    % Workflow config for PTMdecoder top-level orchestration.

    properties
        param_file_path

        % Main peptide workflow mode:
        %   'none' | 'msms_peptide' | 'peptide_requant' | 'peptide_only'
        msms_workflow_mode

        % Independent workflows
        site_level_on
        merge_to_pair_level_on
        merge_pairs_level_on

        % Ordered workflow stages
        execution_order

        % Stage-specific configs
        site_level_config
        merge_to_pair_config
        merge_pairs_config

        % Temporary bridge for legacy modules
        legacy_task_param
    end

    methods
        function obj = CPTMdecoderWorkflowConfig(cfg)
            % Construct from struct config.
            if nargin < 1
                cfg = struct();
            end

            cfg = CPTMdecoderWorkflowConfig.finalize(cfg);

            obj.param_file_path = cfg.param_file_path;
            obj.msms_workflow_mode = cfg.msms_workflow_mode;
            obj.site_level_on = cfg.site_level_on;
            obj.merge_to_pair_level_on = cfg.merge_to_pair_level_on;
            obj.merge_pairs_level_on = cfg.merge_pairs_level_on;
            obj.execution_order = cfg.execution_order;
            obj.site_level_config = cfg.site_level_config;
            obj.merge_to_pair_config = cfg.merge_to_pair_config;
            obj.merge_pairs_config = cfg.merge_pairs_config;
            obj.legacy_task_param = cfg.legacy_task_param;
        end
    end

    methods (Static)
        function obj = fromParamFile(param_file_path)
            % Build workflow config from parameter file.
            taskParam = CTaskParam(param_file_path);
            obj = CPTMdecoderWorkflowConfig.fromTaskParam(taskParam, param_file_path);
        end

        function obj = fromTaskParam(taskParam, param_file_path)
            % Build workflow config from CTaskParam.
            if nargin < 2
                param_file_path = '';
            end

            cfg = struct();
            cfg.param_file_path = param_file_path;
            cfg.legacy_task_param = taskParam;

            cfg.msms_workflow_mode = CPTMdecoderWorkflowConfig.resolveMsmsWorkflowMode(taskParam);
            cfg.site_level_on = logical(taskParam.m_site_level_on);
            cfg.merge_to_pair_level_on = logical(taskParam.m_merge_to_pair_level_on);
            cfg.merge_pairs_level_on = logical(taskParam.m_merge_pairs_level_on);

            cfg.execution_order = { ...
                'msms_workflow', ...
                'site_level', ...
                'merge_to_pair_level', ...
                'merge_pairs_level' ...
            };

            if cfg.site_level_on
                cfg.site_level_config = CSiteLevelPipelineConfig.fromTaskParam(taskParam);
            else
                cfg.site_level_config = [];
            end

            if cfg.merge_to_pair_level_on
                pair_configs = cell(taskParam.m_left_right_out_num, 1);
                for idx_pairs = 1:taskParam.m_left_right_out_num
                    current_pair = taskParam.m_left_right_out(idx_pairs, :);
                    pair_configs{idx_pairs} = CMergeEachPairConfig.fromPairRow( ...
                        current_pair, taskParam.m_left_name, taskParam.m_right_name, ...
                        taskParam.m_ignore_strings_pair_level);
                end
                cfg.merge_to_pair_config = pair_configs;
            else
                cfg.merge_to_pair_config = {};
            end

            if cfg.merge_pairs_level_on
                cfg.merge_pairs_config = CMergePairsConfig.fromTaskParam(taskParam);
            else
                cfg.merge_pairs_config = [];
            end

            obj = CPTMdecoderWorkflowConfig(cfg);
        end
    end

    methods (Static, Access = private)
        function mode = resolveMsmsWorkflowMode(taskParam)
            if taskParam.m_msms_peptide_level_on
                mode = 'msms_peptide';
            elseif taskParam.m_peptide_requant_on
                mode = 'peptide_requant';
            elseif taskParam.m_peptide_only_on
                mode = 'peptide_only';
            else
                mode = 'none';
            end
        end

        function cfg = finalize(cfg)
            if ~isfield(cfg, 'param_file_path') || isempty(cfg.param_file_path)
                cfg.param_file_path = '';
            end
            if ~isfield(cfg, 'msms_workflow_mode') || isempty(cfg.msms_workflow_mode)
                cfg.msms_workflow_mode = 'none';
            end
            if ~isfield(cfg, 'site_level_on') || isempty(cfg.site_level_on)
                cfg.site_level_on = false;
            end
            if ~isfield(cfg, 'merge_to_pair_level_on') || isempty(cfg.merge_to_pair_level_on)
                cfg.merge_to_pair_level_on = false;
            end
            if ~isfield(cfg, 'merge_pairs_level_on') || isempty(cfg.merge_pairs_level_on)
                cfg.merge_pairs_level_on = false;
            end
            if ~isfield(cfg, 'execution_order') || isempty(cfg.execution_order)
                cfg.execution_order = {'msms_workflow', 'site_level', 'merge_to_pair_level', 'merge_pairs_level'};
            end
            if ~isfield(cfg, 'site_level_config') || isempty(cfg.site_level_config)
                cfg.site_level_config = [];
            end
            if ~isfield(cfg, 'merge_to_pair_config') || isempty(cfg.merge_to_pair_config)
                cfg.merge_to_pair_config = {};
            end
            if ~isfield(cfg, 'merge_pairs_config') || isempty(cfg.merge_pairs_config)
                cfg.merge_pairs_config = [];
            end
            if ~isfield(cfg, 'legacy_task_param')
                cfg.legacy_task_param = [];
            end

            valid_modes = {'none', 'msms_peptide', 'peptide_requant', 'peptide_only'};
            if ~any(strcmp(cfg.msms_workflow_mode, valid_modes))
                error('CPTMdecoderWorkflowConfig:InvalidMsmsMode', ...
                    'msms_workflow_mode must be one of: %s', strjoin(valid_modes, ', '));
            end

            if isempty(cfg.legacy_task_param)
                error('CPTMdecoderWorkflowConfig:MissingLegacyTaskParam', ...
                    'legacy_task_param must be provided in current migration phase.');
            end

            if cfg.site_level_on && isempty(cfg.site_level_config)
                error('CPTMdecoderWorkflowConfig:MissingSiteConfig', ...
                    'site_level_config must be provided when site_level_on is true.');
            end

            if cfg.merge_to_pair_level_on && isempty(cfg.merge_to_pair_config)
                error('CPTMdecoderWorkflowConfig:MissingMergeToPairConfig', ...
                    'merge_to_pair_config must be provided when merge_to_pair_level_on is true.');
            end

            if cfg.merge_pairs_level_on && isempty(cfg.merge_pairs_config)
                error('CPTMdecoderWorkflowConfig:MissingMergePairsConfig', ...
                    'merge_pairs_config must be provided when merge_pairs_level_on is true.');
            end
        end
    end
end
