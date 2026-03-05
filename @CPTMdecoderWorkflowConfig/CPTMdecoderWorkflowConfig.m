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

    end

    properties
        param_file_path

        % Ordered workflow stages (execution order is defined by stage order)
        % Each stage is a struct with fields:
        %   - name   : stage name
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

        function flags = resolveMsmsWorkflowFlagsFromMap(task_param_map)
            % Resolve base MSMS workflow flags based on the parameter map.
            %   task_param_map (containers.Map)
            %       parameter key-value map
            msms_pep_on = CPTMdecoderWorkflowConfig.getFlag(task_param_map, CPTMdecoderWorkflowParamKeys.PARAM_MSMS_PEPTIDE_LEVEL_ON);
            msms_only_on = CPTMdecoderWorkflowConfig.getFlag(task_param_map, CPTMdecoderWorkflowParamKeys.PARAM_MSMS_ONLY_ON);
            peptide_only_on = CPTMdecoderWorkflowConfig.getFlag(task_param_map, CPTMdecoderWorkflowParamKeys.PARAM_PEPTIDE_ONLY_ON);

            if (msms_pep_on && msms_only_on) || (msms_pep_on && peptide_only_on) || (msms_only_on && peptide_only_on)
                error('CPTMdecoderWorkflowConfig:InvalidMsmsWorkflowFlags', ...
                    ['msms_peptide_level_on, msms_only_on, peptide_only_on are mutually exclusive ', ...
                    '(msms_peptide_level_on is a legacy combined flag).']);
            end

            flags = struct();
            if msms_pep_on
                warning('CPTMdecoderWorkflowConfig:DeprecatedMsmsPeptideLevelFlag', ...
                    ['msms_peptide_level_on is deprecated. Please use msms_only_on=1 and peptide_only_on=1 ', ...
                    'in separate runs when needed.']);
                flags.run_msms_level = true;
                flags.run_peptide_quant = true;
            else
                flags.run_msms_level = msms_only_on;
                flags.run_peptide_quant = peptide_only_on;
            end
        end

        function stage = makeStage(name, cfg_obj, enabled)
            % Build a stage struct.
            %   name (1 x 1 char/string)
            %       stage name
            %   cfg_obj (struct)
            %       config object
            %   enabled (logical)
            %       whether the stage is enabled
            stage = struct();
            stage.name = name;
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
                required_fields = {'name', 'config', 'enabled'};
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
            end
        end
    end
end
