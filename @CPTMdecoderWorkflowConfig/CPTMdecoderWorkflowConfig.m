classdef CPTMdecoderWorkflowConfig
    % Workflow config for PTMdecoder top-level orchestration.

    properties (Constant)
        STAGE_MSMS_LEVEL = 'msms_level'
        STAGE_PEPTIDE_QUANT = 'peptide_quant'
        STAGE_PEPTIDE_REQUANT = 'peptide_requant'
        STAGE_PEPTIDE_ALIGN_REQUANT = 'peptide_align_requant'
        STAGE_NORM_PEPTIDE_QUANT = 'norm_peptide_quant'
        STAGE_NORM_PEPTIDE_REQUANT = 'norm_peptide_requant'
        STAGE_SITE_LEVEL = 'site_level'
        STAGE_MERGE_TO_PAIR_LEVEL = 'merge_to_pair_level'
        STAGE_MERGE_PAIRS_LEVEL = 'merge_pairs_level'

        ACTION_MSMS_PEPTIDE = 'msms_peptide'
        ACTION_MSMS_ONLY = 'msms_only'
        ACTION_PEPTIDE_REQUANT = 'peptide_requant'
        ACTION_PEPTIDE_ALIGN_REQUANT = 'peptide_align_requant'
        ACTION_PEPTIDE_ONLY = 'peptide_only'
        ACTION_NORM_PEPTIDE_QUANT = 'norm_peptide_quant'
        ACTION_NORM_PEPTIDE_REQUANT = 'norm_peptide_requant'
        ACTION_SITE_SUMMARY = 'site_summary'
        ACTION_MERGE_EACH_PAIR = 'merge_each_pair'
        ACTION_MERGE_PAIRS = 'merge_pairs'

    end

    properties
        param_file_path

        % Ordered workflow stages (execution order is defined by stage order)
        % Each stage is a struct with fields:
        %   - name   : stage name
        %   - action : compatibility field for input intent trace; not used for runtime dispatch
        %   - config : stage config object
        %   - enabled: logical
        stages
    end

    methods
        function obj = CPTMdecoderWorkflowConfig(cfg)
            % Construct from struct config.
            if nargin < 1
                cfg = struct();
            end

            cfg = CPTMdecoderWorkflowConfig.finalize(cfg);

            obj.param_file_path = cfg.param_file_path;
            obj.stages = cfg.stages;
        end
    end

    methods (Static)
        % Build workflow config directly from parameter file.
        obj = fromParamFile(param_file_path)
    end

    methods (Static, Access = private)
        function flag = getFlag(task_param_map, key_name)
            % Get a boolean flag from the parameter map.
            %   task_param_map (containers.Map)
            %       parameter key-value map
            %   key_name (1 x 1 char/string)
            %       parameter key name
            flag = task_param_map.isKey(key_name) && strcmp(strtrim(task_param_map(key_name)), '1');
        end

        function action = resolveMsmsWorkflowActionFromMap(task_param_map)
            % Resolve the MSMS workflow action based on flags in the parameter map.
            %   task_param_map (containers.Map)
            %       parameter key-value map
            msms_pep_on = CPTMdecoderWorkflowConfig.getFlag(task_param_map, CPTMdecoderWorkflowParamKeys.PARAM_MSMS_PEPTIDE_LEVEL_ON);
            msms_only_on = CPTMdecoderWorkflowConfig.getFlag(task_param_map, CPTMdecoderWorkflowParamKeys.PARAM_MSMS_ONLY_ON);
            requant_on = CPTMdecoderWorkflowConfig.getFlag(task_param_map, CPTMdecoderWorkflowParamKeys.PARAM_PEPTIDE_REQUANT_ON);
            align_requant_on = CPTMdecoderWorkflowConfig.getFlag(task_param_map, CPTMdecoderWorkflowParamKeys.PARAM_PEPTIDE_ALIGN_REQUANT_ON);
            peptide_only_on = CPTMdecoderWorkflowConfig.getFlag(task_param_map, CPTMdecoderWorkflowParamKeys.PARAM_PEPTIDE_ONLY_ON);

            if msms_pep_on + msms_only_on + requant_on + align_requant_on + peptide_only_on > 1
                error('CPTMdecoderWorkflowConfig:InvalidMsmsWorkflowAction', ...
                    ['msms_peptide_level_on, msms_only_on, peptide_requant_on, ', ...
                    'peptide_align_requant_on, peptide_only_on cannot be enabled simultaneously.']);
            end

            if msms_pep_on
                action = CPTMdecoderWorkflowConfig.ACTION_MSMS_PEPTIDE;
            elseif msms_only_on
                action = CPTMdecoderWorkflowConfig.ACTION_MSMS_ONLY;
            elseif requant_on
                action = CPTMdecoderWorkflowConfig.ACTION_PEPTIDE_REQUANT;
            elseif align_requant_on
                action = CPTMdecoderWorkflowConfig.ACTION_PEPTIDE_ALIGN_REQUANT;
            elseif peptide_only_on
                action = CPTMdecoderWorkflowConfig.ACTION_PEPTIDE_ONLY;
            else
                action = '';
            end
        end

        function stage = makeStage(name, action, cfg_obj, enabled)
            % Build a stage struct.
            %   name (1 x 1 char/string)
            %       stage name
            %   action (1 x 1 char/string)
            %       stage action
            %   cfg_obj (struct)
            %       config object
            %   enabled (logical)
            %       whether the stage is enabled
            % Build a stage struct.
            %   name (1 x 1 char/string)
            %       stage name
            %   action (1 x 1 char/string)
            %       stage action
            %   cfg_obj (struct)
            %       config object
            %   enabled (logical)
            %       whether the stage is enabled
            % Build a stage struct.
            % Note: action is retained as a compatibility field for input intent trace,
            % while runtime dispatch is based on stage.name.
            stage = struct();
            stage.name = name;
            stage.action = action;
            stage.config = cfg_obj;
            stage.enabled = logical(enabled);
        end

        function names = allStageNames()
            % Get all valid stage names.
            names = { ...
                CPTMdecoderWorkflowConfig.STAGE_MSMS_LEVEL, ...
                CPTMdecoderWorkflowConfig.STAGE_PEPTIDE_QUANT, ...
                CPTMdecoderWorkflowConfig.STAGE_PEPTIDE_REQUANT, ...
                CPTMdecoderWorkflowConfig.STAGE_PEPTIDE_ALIGN_REQUANT, ...
                CPTMdecoderWorkflowConfig.STAGE_NORM_PEPTIDE_QUANT, ...
                CPTMdecoderWorkflowConfig.STAGE_NORM_PEPTIDE_REQUANT, ...
                CPTMdecoderWorkflowConfig.STAGE_SITE_LEVEL, ...
                CPTMdecoderWorkflowConfig.STAGE_MERGE_TO_PAIR_LEVEL, ...
                CPTMdecoderWorkflowConfig.STAGE_MERGE_PAIRS_LEVEL ...
            };
        end

        function actions = stageActions(stage_name)
            % Get valid actions for a given stage name.
            %   stage_name (1 x 1 char/string)
            %       stage name
            switch stage_name
                case CPTMdecoderWorkflowConfig.STAGE_MSMS_LEVEL
                    actions = {CPTMdecoderWorkflowConfig.ACTION_MSMS_ONLY};
                case CPTMdecoderWorkflowConfig.STAGE_PEPTIDE_QUANT
                    actions = {CPTMdecoderWorkflowConfig.ACTION_PEPTIDE_ONLY};
                case CPTMdecoderWorkflowConfig.STAGE_PEPTIDE_REQUANT
                    actions = {CPTMdecoderWorkflowConfig.ACTION_PEPTIDE_REQUANT};
                case CPTMdecoderWorkflowConfig.STAGE_PEPTIDE_ALIGN_REQUANT
                    actions = {CPTMdecoderWorkflowConfig.ACTION_PEPTIDE_ALIGN_REQUANT};
                case CPTMdecoderWorkflowConfig.STAGE_NORM_PEPTIDE_QUANT
                    actions = {CPTMdecoderWorkflowConfig.ACTION_NORM_PEPTIDE_QUANT};
                case CPTMdecoderWorkflowConfig.STAGE_NORM_PEPTIDE_REQUANT
                    actions = {CPTMdecoderWorkflowConfig.ACTION_NORM_PEPTIDE_REQUANT};
                case CPTMdecoderWorkflowConfig.STAGE_SITE_LEVEL
                    actions = {CPTMdecoderWorkflowConfig.ACTION_SITE_SUMMARY};
                case CPTMdecoderWorkflowConfig.STAGE_MERGE_TO_PAIR_LEVEL
                    actions = {CPTMdecoderWorkflowConfig.ACTION_MERGE_EACH_PAIR};
                case CPTMdecoderWorkflowConfig.STAGE_MERGE_PAIRS_LEVEL
                    actions = {CPTMdecoderWorkflowConfig.ACTION_MERGE_PAIRS};
                otherwise
                    actions = {};
            end
        end

        function cfg = finalize(cfg)
            % Finalize and validate the config struct.
            %   cfg (struct)
            %       config struct to finalize
            if ~isfield(cfg, 'param_file_path') || isempty(cfg.param_file_path)
                cfg.param_file_path = '';
            end
            if ~isfield(cfg, 'stages') || isempty(cfg.stages)
                cfg.stages = {};
            end

            if ~iscell(cfg.stages)
                error('CPTMdecoderWorkflowConfig:InvalidStages', ...
                    'stages must be a cell array of stage structs.');
            end

            for i_stage = 1:numel(cfg.stages)
                stage = cfg.stages{i_stage};
                if ~isstruct(stage)
                    error('CPTMdecoderWorkflowConfig:InvalidStageType', ...
                        'Each stage must be a struct.');
                end
                required_fields = {'name', 'action', 'config', 'enabled'};
                for i_field = 1:numel(required_fields)
                    if ~isfield(stage, required_fields{i_field})
                        error('CPTMdecoderWorkflowConfig:InvalidStageField', ...
                            'Each stage must include field ''%s''.', required_fields{i_field});
                    end
                end
                if ~islogical(stage.enabled)
                    error('CPTMdecoderWorkflowConfig:InvalidStageEnabled', ...
                        'stage.enabled must be logical.');
                end
                if ~isscalar(stage.enabled)
                    error('CPTMdecoderWorkflowConfig:InvalidStageEnabledSize', ...
                        'stage.enabled must be a logical scalar.');
                end
                valid_names = CPTMdecoderWorkflowConfig.allStageNames();
                if ~any(strcmp(stage.name, valid_names))
                    error('CPTMdecoderWorkflowConfig:InvalidStageName', ...
                        'Unsupported stage.name: %s', stage.name);
                end
                valid_actions = CPTMdecoderWorkflowConfig.stageActions(stage.name);
                if ~any(strcmp(stage.action, valid_actions))
                    error('CPTMdecoderWorkflowConfig:InvalidStageAction', ...
                        ['Unsupported compatibility action ''%s'' for stage ''%s'' ', ...
                        '(runtime dispatch is based on stage.name).'], stage.action, stage.name);
                end
            end
        end
    end
end
