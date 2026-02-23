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

            cfg.site_level_config = struct( ...
                'input_path', taskParam.m_pep_level_file_path, ...
                'output_intere_path', taskParam.m_output_intere_path, ...
                'output_unintere_path', taskParam.m_output_unintere_path, ...
                'protein_name_abbr', taskParam.m_protein_name_abbr, ...
                'mod_name_abbr', taskParam.m_mod_name_abbr, ...
                'ignore_strings', {taskParam.m_ignore_strings_site_level}, ...
                'column_idxs', struct('icol_seq', 2, 'icol_auc', 8));

            cfg.merge_to_pair_config = struct( ...
                'pairs', {taskParam.m_left_right_out}, ...
                'left_name', taskParam.m_left_name, ...
                'right_name', taskParam.m_right_name, ...
                'ignore_strings', {taskParam.m_ignore_strings_pair_level});

            cfg.merge_pairs_config = struct( ...
                'result_paths', {taskParam.m_pair}, ...
                'output_path', taskParam.m_final_output_path, ...
                'group_titles', {taskParam.m_left_right_name});

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
                cfg.site_level_config = struct();
            end
            if ~isfield(cfg, 'merge_to_pair_config') || isempty(cfg.merge_to_pair_config)
                cfg.merge_to_pair_config = struct('pairs', {cell(0, 3)}, 'left_name', '', 'right_name', '', 'ignore_strings', {{}});
            end
            if ~isfield(cfg, 'merge_pairs_config') || isempty(cfg.merge_pairs_config)
                cfg.merge_pairs_config = struct('result_paths', {{}}, 'output_path', '', 'group_titles', {{}});
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
        end
    end
end
