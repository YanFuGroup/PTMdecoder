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
        msms_peptide_config
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
            obj.msms_peptide_config = cfg.msms_peptide_config;
        end
    end

    methods (Static)
        function obj = fromParamFile(param_file_path)
            % Build workflow config directly from parameter file.
            task_param_map = CPTMdecoderWorkflowConfig.parseParamFileToMap(param_file_path);

            cfg = struct();
            cfg.param_file_path = param_file_path;

            cfg.msms_workflow_mode = CPTMdecoderWorkflowConfig.resolveMsmsWorkflowModeFromMap(task_param_map);
            cfg.site_level_on = CPTMdecoderWorkflowConfig.getFlag(task_param_map, 'site_level_on');
            cfg.merge_to_pair_level_on = CPTMdecoderWorkflowConfig.getFlag(task_param_map, 'merge_to_pair_level_on');
            cfg.merge_pairs_level_on = CPTMdecoderWorkflowConfig.getFlag(task_param_map, 'merge_pairs_level_on');

            cfg.execution_order = { ...
                'msms_workflow', ...
                'site_level', ...
                'merge_to_pair_level', ...
                'merge_pairs_level' ...
            };

            if ~strcmp(cfg.msms_workflow_mode, 'none')
                msms_cfg_struct = CPTMdecoderWorkflowConfig.buildMsmsPepConfigStructFromMap(task_param_map, cfg.msms_workflow_mode);
                cfg.msms_peptide_config = CMSMSPepDeconvConfig(msms_cfg_struct);
            else
                cfg.msms_peptide_config = [];
            end

            if cfg.site_level_on
                site_cfg_struct = CPTMdecoderWorkflowConfig.buildSiteConfigStructFromMap(task_param_map);
                cfg.site_level_config = CSiteLevelPipelineConfig(site_cfg_struct);
            else
                cfg.site_level_config = [];
            end

            if cfg.merge_to_pair_level_on
                cfg.merge_to_pair_config = CPTMdecoderWorkflowConfig.buildMergeToPairConfigsFromMap(task_param_map);
            else
                cfg.merge_to_pair_config = {};
            end

            if cfg.merge_pairs_level_on
                merge_pairs_cfg_struct = CPTMdecoderWorkflowConfig.buildMergePairsConfigStructFromMap(task_param_map);
                cfg.merge_pairs_config = CMergePairsConfig(merge_pairs_cfg_struct);
            else
                cfg.merge_pairs_config = [];
            end

            obj = CPTMdecoderWorkflowConfig(cfg);
        end

    end

    methods (Static, Access = private)
        function task_param_map = parseParamFileToMap(param_file_path)
            fid = fopen(param_file_path, 'r');
            if fid <= 0
                error('CPTMdecoderWorkflowConfig:OpenParamFileFailed', ...
                    'Failed to open parameter file: %s', param_file_path);
            end

            task_param_map = containers.Map();
            line_num = 0;
            while ~feof(fid)
                line_num = line_num + 1;
                str_line = fgetl(fid);
                if ~ischar(str_line)
                    continue;
                end

                str_line = CPTMdecoderWorkflowConfig.removeComments(str_line);
                str_line = strtrim(str_line);
                if isempty(str_line)
                    continue;
                end

                str_seg = split(str_line, '=');
                if length(str_seg) ~= 2
                    error('CPTMdecoderWorkflowConfig:InvalidParamFormat', ...
                        'Unexpected parameter format in line %d: %s', line_num, str_line);
                end

                m_key = strtrim(str_seg(1));
                m_val = strtrim(str_seg(2));
                task_param_map(m_key{1}) = m_val{1};
            end
            fclose(fid);
        end

        function str_line = removeComments(str_line)
            idx = strfind(str_line, '#');
            if ~isempty(idx)
                str_line(idx(1):end) = [];
            end
        end

        function flag = getFlag(task_param_map, key_name)
            flag = task_param_map.isKey(key_name) && strcmp(strtrim(task_param_map(key_name)), '1');
        end

        function value = getRequired(task_param_map, key_name, explain)
            if ~task_param_map.isKey(key_name)
                error('CPTMdecoderWorkflowConfig:MissingRequiredParam', ...
                    'Required param ''%s'' is missing (%s).', key_name, explain);
            end
            value = task_param_map(key_name);
        end

        function value = getOptional(task_param_map, key_name, default_value)
            if task_param_map.isKey(key_name)
                value = task_param_map(key_name);
            else
                value = default_value;
            end
        end

        function mode = resolveMsmsWorkflowModeFromMap(task_param_map)
            msms_on = CPTMdecoderWorkflowConfig.getFlag(task_param_map, 'msms_peptide_level_on');
            requant_on = CPTMdecoderWorkflowConfig.getFlag(task_param_map, 'peptide_requant_on');
            peptide_only_on = CPTMdecoderWorkflowConfig.getFlag(task_param_map, 'peptide_only_on');

            if msms_on + requant_on + peptide_only_on > 1
                error('CPTMdecoderWorkflowConfig:InvalidMsmsWorkflowMode', ...
                    'msms_peptide_level_on, peptide_requant_on, peptide_only_on cannot be enabled simultaneously.');
            end

            if msms_on
                mode = 'msms_peptide';
            elseif requant_on
                mode = 'peptide_requant';
            elseif peptide_only_on
                mode = 'peptide_only';
            else
                mode = 'none';
            end
        end

        function cfg = buildMsmsPepConfigStructFromMap(task_param_map, mode)
            cfg = struct();
            cfg.mod_file_path = CPTMdecoderWorkflowConfig.getRequired(task_param_map, 'mod_file_path', 'modification file path');
            cfg.fixed_mod = CPTMdecoderWorkflowConfig.getRequired(task_param_map, 'fixed_mod', 'fixed modification list');
            cfg.variable_mod = CPTMdecoderWorkflowConfig.getRequired(task_param_map, 'variable_mod', 'variable modification list');
            cfg.spec_dir_path = CPTMdecoderWorkflowConfig.getRequired(task_param_map, 'spec_dir_path', 'spectrum directory path');

            ms1_tol_value = str2double(CPTMdecoderWorkflowConfig.getRequired(task_param_map, 'ms1_tolerance_value', 'MS1 tolerance value'));
            ms1_tol_type = strtrim(CPTMdecoderWorkflowConfig.getRequired(task_param_map, 'ms1_tolerance_type', 'MS1 tolerance type'));
            cfg.ms1_tolerance = struct('value', ms1_tol_value, 'isppm', strcmpi(ms1_tol_type, 'PPM'));

            cfg.ms2_tolerance = str2double(CPTMdecoderWorkflowConfig.getRequired(task_param_map, 'ms2_tolerance', 'MS2 tolerance value'));
            cfg.alpha = str2double(CPTMdecoderWorkflowConfig.getRequired(task_param_map, 'alpha', 'noise filtering threshold'));
            cfg.fasta_file_path = CPTMdecoderWorkflowConfig.getRequired(task_param_map, 'fasta_file_path', 'fasta file path');
            cfg.regular_express = CPTMdecoderWorkflowConfig.getRequired(task_param_map, 'regular_express', 'protein name regex');
            cfg.filtered_res_file_path = CPTMdecoderWorkflowConfig.getOptional(task_param_map, 'filtered_res_file_path', '');
            cfg.model = str2double(CPTMdecoderWorkflowConfig.getRequired(task_param_map, 'model', 'algorithm model'));
            cfg.method = str2double(CPTMdecoderWorkflowConfig.getRequired(task_param_map, 'method', 'algorithm method'));

            lambda_value = CPTMdecoderWorkflowConfig.getOptional(task_param_map, 'lambda', '0');
            cfg.lambda = str2double(lambda_value);
            if isequal(strtrim(CPTMdecoderWorkflowConfig.getOptional(task_param_map, 'method', '')), '2') && ~task_param_map.isKey('lambda')
                error('CPTMdecoderWorkflowConfig:MissingLambda', ...
                    'The lasso parameter ''lambda'' is required when method=2.');
            end

            cfg.result_filter_threshold = str2double(CPTMdecoderWorkflowConfig.getRequired(task_param_map, 'result_filter_threshold', 'result filter threshold'));
            cfg.enzyme_name = CPTMdecoderWorkflowConfig.getRequired(task_param_map, 'enzyme_name', 'enzyme name');
            cfg.enzyme_limits = str2num(CPTMdecoderWorkflowConfig.getRequired(task_param_map, 'enzyme_limit_C_term_possible_mod', 'enzyme C-term possible modifications')); %#ok<ST2NM>
            cfg.output_dir_path = CPTMdecoderWorkflowConfig.getRequired(task_param_map, 'output_dir_path', 'output directory');

            cfg.min_MSMS_num = str2double(CPTMdecoderWorkflowConfig.getOptional(task_param_map, 'min_MSMS_num', '1'));
            cfg.ion_types = [1,2];
            cfg.case_penalty_intens = 'intens_sum';
            cfg.grid_penalty_intens = 'intens_sum';
            cfg.case_OLS_intens_weight = 'none';

            if strcmp(mode, 'msms_peptide')
                cfg.pep_spec_file_path = CPTMdecoderWorkflowConfig.getRequired(task_param_map, 'pep_spec_file_path', 'peptide-spectrum list path');
                cfg.checked_peptides_res_path = [];
                cfg.msms_res_path = [];
            elseif strcmp(mode, 'peptide_requant')
                cfg.pep_spec_file_path = '';
                cfg.checked_peptides_res_path = CPTMdecoderWorkflowConfig.getRequired(task_param_map, 'checked_peptides_res_path', 'checked peptide result path');
                cfg.msms_res_path = CPTMdecoderWorkflowConfig.getRequired(task_param_map, 'msms_res_path', 'MSMS result path');
            else
                cfg.pep_spec_file_path = '';
                cfg.checked_peptides_res_path = [];
                cfg.msms_res_path = CPTMdecoderWorkflowConfig.getRequired(task_param_map, 'msms_res_path', 'MSMS result path');
            end
        end

        function cfg = buildSiteConfigStructFromMap(task_param_map)
            cfg = struct();

            output_dir_path = CPTMdecoderWorkflowConfig.getOptional(task_param_map, 'output_dir_path', '');
            pep_default = fullfile(output_dir_path, 'report_peptide_all.txt');
            intere_default = fullfile(output_dir_path, 'report_site.txt');
            unintere_default = fullfile(output_dir_path, 'report_peptide_uninterested.txt');

            cfg.input_path = CPTMdecoderWorkflowConfig.getOptional(task_param_map, 'pep_level_file_path', pep_default);
            if isempty(cfg.input_path)
                cfg.input_path = pep_default;
            end

            cfg.output_intere_path = CPTMdecoderWorkflowConfig.getOptional(task_param_map, 'output_intere_path', intere_default);
            if isempty(cfg.output_intere_path)
                cfg.output_intere_path = intere_default;
            end

            cfg.output_unintere_path = CPTMdecoderWorkflowConfig.getOptional(task_param_map, 'output_unintere_path', unintere_default);
            if isempty(cfg.output_unintere_path)
                cfg.output_unintere_path = unintere_default;
            end

            protein_name_abbr_num = str2double(CPTMdecoderWorkflowConfig.getRequired(task_param_map, 'protein_name_abbr_num', 'number of protein abbreviation mappings'));
            protein_name_abbr = containers.Map;
            for idx = 1:protein_name_abbr_num
                key_name = ['protein_name_abbr_', num2str(idx)];
                pair_str = CPTMdecoderWorkflowConfig.getRequired(task_param_map, key_name, 'protein abbreviation pair');
                split_str = strsplit(pair_str, '>');
                protein_name_abbr(strtrim(split_str{1})) = strtrim(split_str{2});
            end
            cfg.protein_name_abbr = protein_name_abbr;

            mod_name_abbr_num = str2double(CPTMdecoderWorkflowConfig.getRequired(task_param_map, 'mod_name_abbr_num', 'number of modification abbreviation mappings'));
            mod_name_abbr = containers.Map;
            for idx = 1:mod_name_abbr_num
                key_name = ['mod_name_abbr_', num2str(idx)];
                pair_str = CPTMdecoderWorkflowConfig.getRequired(task_param_map, key_name, 'modification abbreviation pair');
                split_str = strsplit(pair_str, '>');
                mod_name_abbr(strtrim(split_str{1})) = strtrim(split_str{2});
            end
            cfg.mod_name_abbr = mod_name_abbr;

            ignore_str = CPTMdecoderWorkflowConfig.getRequired(task_param_map, 'ignore_strings_site_level', 'ignore strings for site-level summary');
            cfg.ignore_strings = CPTMdecoderWorkflowConfig.parseIgnoreStrings(ignore_str);
            cfg.column_idxs = struct('icol_seq', 2, 'icol_auc', 8);
        end

        function pair_configs = buildMergeToPairConfigsFromMap(task_param_map)
            left_right_out_num = str2double(CPTMdecoderWorkflowConfig.getRequired(task_param_map, 'left_right_out_num', 'number of pairwise comparisons'));
            left_name = CPTMdecoderWorkflowConfig.getRequired(task_param_map, 'left_name', 'left group name');
            right_name = CPTMdecoderWorkflowConfig.getRequired(task_param_map, 'right_name', 'right group name');
            ignore_str = CPTMdecoderWorkflowConfig.getRequired(task_param_map, 'ignore_strings_pair_level', 'ignore strings for pair-level merge');
            ignore_strings = CPTMdecoderWorkflowConfig.parseIgnoreStrings(ignore_str);

            pair_configs = cell(left_right_out_num, 1);
            for idx = 1:left_right_out_num
                key_name = ['left_right_out_', num2str(idx)];
                left_right_out_str = CPTMdecoderWorkflowConfig.getRequired(task_param_map, key_name, 'pairwise input/output mapping');
                split_str = strsplit(left_right_out_str, {'|', '>'});
                pair_row = strtrim(split_str);
                pair_configs{idx} = CMergeEachPairConfig.fromPairRow(pair_row, left_name, right_name, ignore_strings);
            end
        end

        function cfg = buildMergePairsConfigStructFromMap(task_param_map)
            pair_num = str2double(CPTMdecoderWorkflowConfig.getRequired(task_param_map, 'pair_num', 'number of pairs to merge'));
            result_paths = cell(pair_num, 1);
            group_titles = cell(pair_num, 2);
            for idx = 1:pair_num
                pair_key = ['pair_', num2str(idx)];
                lr_key = ['left_right_name_', num2str(idx)];
                result_paths{idx} = CPTMdecoderWorkflowConfig.getRequired(task_param_map, pair_key, 'pair-level result path');

                lr_str = CPTMdecoderWorkflowConfig.getRequired(task_param_map, lr_key, 'left-right group names');
                split_lr = strsplit(lr_str, '|');
                group_titles{idx, 1} = strtrim(split_lr{1});
                group_titles{idx, 2} = strtrim(split_lr{2});
            end

            cfg = struct();
            cfg.result_paths = result_paths;
            cfg.output_path = CPTMdecoderWorkflowConfig.getRequired(task_param_map, 'final_output_path', 'final merged output path');
            cfg.group_titles = group_titles;
        end

        function ignore_strings = parseIgnoreStrings(ignore_strings_str)
            ignore_strings = {};
            ignore_strings_seg = strsplit(ignore_strings_str, ';');
            ignore_strings_seg = cellfun(@(s) regexp(s, '"([^"]*)"', 'tokens'), ...
                ignore_strings_seg, 'UniformOutput', false);
            for idx_is = 1:length(ignore_strings_seg)
                if isempty(ignore_strings_seg{idx_is})
                    continue;
                end
                ignore_strings = [ignore_strings, ignore_strings_seg{idx_is}{1}]; %#ok<AGROW>
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
            if ~isfield(cfg, 'msms_peptide_config')
                cfg.msms_peptide_config = [];
            end

            valid_modes = {'none', 'msms_peptide', 'peptide_requant', 'peptide_only'};
            if ~any(strcmp(cfg.msms_workflow_mode, valid_modes))
                error('CPTMdecoderWorkflowConfig:InvalidMsmsMode', ...
                    'msms_workflow_mode must be one of: %s', strjoin(valid_modes, ', '));
            end

            if ~strcmp(cfg.msms_workflow_mode, 'none') && isempty(cfg.msms_peptide_config)
                error('CPTMdecoderWorkflowConfig:MissingMsmsConfig', ...
                    'msms_peptide_config must be provided when msms_workflow_mode is not none.');
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
