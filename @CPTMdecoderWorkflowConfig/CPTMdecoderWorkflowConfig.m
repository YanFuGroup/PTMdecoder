classdef CPTMdecoderWorkflowConfig
    % Workflow config for PTMdecoder top-level orchestration.

    properties (Constant)
        STAGE_MSMS_LEVEL = 'msms_level'
        STAGE_PEPTIDE_QUANT = 'peptide_quant'
        STAGE_PEPTIDE_REQUANT = 'peptide_requant'
        STAGE_NORM_PEPTIDE_QUANT = 'norm_peptide_quant'
        STAGE_NORM_PEPTIDE_REQUANT = 'norm_peptide_requant'
        STAGE_SITE_LEVEL = 'site_level'
        STAGE_MERGE_TO_PAIR_LEVEL = 'merge_to_pair_level'
        STAGE_MERGE_PAIRS_LEVEL = 'merge_pairs_level'

        ACTION_MSMS_PEPTIDE = 'msms_peptide'
        ACTION_MSMS_ONLY = 'msms_only'
        ACTION_PEPTIDE_REQUANT = 'peptide_requant'
        ACTION_PEPTIDE_ONLY = 'peptide_only'
        ACTION_NORM_PEPTIDE_QUANT = 'norm_peptide_quant'
        ACTION_NORM_PEPTIDE_REQUANT = 'norm_peptide_requant'
        ACTION_SITE_SUMMARY = 'site_summary'
        ACTION_MERGE_EACH_PAIR = 'merge_each_pair'
        ACTION_MERGE_PAIRS = 'merge_pairs'

        PARAM_MSMS_PEPTIDE_LEVEL_ON = 'msms_peptide_level_on'
        PARAM_MSMS_ONLY_ON = 'msms_only_on'
        PARAM_PEPTIDE_ONLY_ON = 'peptide_only_on'
        PARAM_PEPTIDE_REQUANT_ON = 'peptide_requant_on'
        PARAM_NORM_PEPTIDE_QUANT_ON = 'norm_peptide_quant_on'
        PARAM_NORM_PEPTIDE_REQUANT_ON = 'norm_peptide_requant_on'
        PARAM_SITE_LEVEL_ON = 'site_level_on'
        PARAM_MERGE_TO_PAIR_LEVEL_ON = 'merge_to_pair_level_on'
        PARAM_MERGE_PAIRS_LEVEL_ON = 'merge_pairs_level_on'

        PARAM_MOD_FILE_PATH = 'mod_file_path'
        PARAM_FIXED_MOD = 'fixed_mod'
        PARAM_VARIABLE_MOD = 'variable_mod'
        PARAM_SPEC_DIR_PATH = 'spec_dir_path'
        PARAM_MS1_TOLERANCE_VALUE = 'ms1_tolerance_value'
        PARAM_MS1_TOLERANCE_TYPE = 'ms1_tolerance_type'
        PARAM_MS2_TOLERANCE = 'ms2_tolerance'
        PARAM_ALPHA = 'alpha'
        PARAM_FASTA_FILE_PATH = 'fasta_file_path'
        PARAM_REGULAR_EXPRESS = 'regular_express'
        PARAM_FILTERED_RES_FILE_PATH = 'filtered_res_file_path'
        PARAM_MODEL = 'model'
        PARAM_METHOD = 'method'
        PARAM_LAMBDA = 'lambda'
        PARAM_RESULT_FILTER_THRESHOLD = 'result_filter_threshold'
        PARAM_ENZYME_NAME = 'enzyme_name'
        PARAM_ENZYME_LIMIT_C_TERM_POSSIBLE_MOD = 'enzyme_limit_C_term_possible_mod'
        PARAM_OUTPUT_DIR_PATH = 'output_dir_path'
        PARAM_MIN_MSMS_NUM = 'min_MSMS_num'
        PARAM_PEP_SPEC_FILE_PATH = 'pep_spec_file_path'
        PARAM_CHECKED_PEPTIDES_RES_PATH = 'checked_peptides_res_path'
        PARAM_MSMS_RES_PATH = 'msms_res_path'
        PARAM_NORM_PROTEIN_PEPTIDE_PAIR_PATH = 'norm_protein_peptide_pair_path'
        PARAM_NORM_OUTPUT_FILE_NAME = 'norm_output_file_name'
        PARAM_NORM_REQUANT_OUTPUT_FILE_NAME = 'norm_requant_output_file_name'

        PARAM_PEP_LEVEL_FILE_PATH = 'pep_level_file_path'
        PARAM_OUTPUT_INTERE_PATH = 'output_intere_path'
        PARAM_OUTPUT_UNINTERE_PATH = 'output_unintere_path'
        PARAM_PROTEIN_NAME_ABBR_NUM = 'protein_name_abbr_num'
        PARAM_MOD_NAME_ABBR_NUM = 'mod_name_abbr_num'
        PARAM_IGNORE_STRINGS_SITE_LEVEL = 'ignore_strings_site_level'

        PARAM_LEFT_RIGHT_OUT_NUM = 'left_right_out_num'
        PARAM_LEFT_NAME = 'left_name'
        PARAM_RIGHT_NAME = 'right_name'
        PARAM_IGNORE_STRINGS_PAIR_LEVEL = 'ignore_strings_pair_level'

        PARAM_PAIR_NUM = 'pair_num'
        PARAM_FINAL_OUTPUT_PATH = 'final_output_path'

        PARAM_PREFIX_PROTEIN_NAME_ABBR = 'protein_name_abbr_'
        PARAM_PREFIX_MOD_NAME_ABBR = 'mod_name_abbr_'
        PARAM_PREFIX_LEFT_RIGHT_OUT = 'left_right_out_'
        PARAM_PREFIX_PAIR = 'pair_'
        PARAM_PREFIX_LEFT_RIGHT_NAME = 'left_right_name_'
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

        function value = getRequired(task_param_map, key_name, explain)
            % Get a required parameter value from the map, error if missing.
            %   task_param_map (containers.Map)
            %       parameter key-value map
            %   key_name (1 x 1 char/string)
            %       parameter key name
            %   explain (1 x 1 char/string)
            %       description for error message
            if ~task_param_map.isKey(key_name)
                error('CPTMdecoderWorkflowConfig:MissingRequiredParam', ...
                    'Required param ''%s'' is missing (%s).', key_name, explain);
            end
            value = task_param_map(key_name);
        end

        function value = getOptional(task_param_map, key_name, default_value)
            % Get an optional parameter value from the map, use default if missing.
            %   task_param_map (containers.Map)
            %       parameter key-value map
            %   key_name (1 x 1 char/string)
            %       parameter key name
            %   default_value (any)
            %       default value if key not present
            if task_param_map.isKey(key_name)
                value = task_param_map(key_name);
            else
                value = default_value;
            end
        end

        function action = resolveMsmsWorkflowActionFromMap(task_param_map)
            % Resolve the MSMS workflow action based on flags in the parameter map.
            %   task_param_map (containers.Map)
            %       parameter key-value map
            msms_pep_on = CPTMdecoderWorkflowConfig.getFlag(task_param_map, CPTMdecoderWorkflowConfig.PARAM_MSMS_PEPTIDE_LEVEL_ON);
            msms_only_on = CPTMdecoderWorkflowConfig.getFlag(task_param_map, CPTMdecoderWorkflowConfig.PARAM_MSMS_ONLY_ON);
            requant_on = CPTMdecoderWorkflowConfig.getFlag(task_param_map, CPTMdecoderWorkflowConfig.PARAM_PEPTIDE_REQUANT_ON);
            peptide_only_on = CPTMdecoderWorkflowConfig.getFlag(task_param_map, CPTMdecoderWorkflowConfig.PARAM_PEPTIDE_ONLY_ON);

            if msms_pep_on + msms_only_on + requant_on + peptide_only_on > 1
                error('CPTMdecoderWorkflowConfig:InvalidMsmsWorkflowAction', ...
                    'msms_peptide_level_on, msms_only_on, peptide_requant_on, peptide_only_on cannot be enabled simultaneously.');
            end

            if msms_pep_on
                action = CPTMdecoderWorkflowConfig.ACTION_MSMS_PEPTIDE;
            elseif msms_only_on
                action = CPTMdecoderWorkflowConfig.ACTION_MSMS_ONLY;
            elseif requant_on
                action = CPTMdecoderWorkflowConfig.ACTION_PEPTIDE_REQUANT;
            elseif peptide_only_on
                action = CPTMdecoderWorkflowConfig.ACTION_PEPTIDE_ONLY;
            else
                action = '';
            end
        end

        function cfg = buildMsmsPepConfigStructFromMap(task_param_map, mode)
            % Build MSMS peptide config struct from parameter map.
            %   task_param_map (containers.Map)
            %       parameter key-value map
            %   mode (1 x 1 char/string)
            %       workflow mode
            cfg = struct();
            cfg.mod_file_path = CPTMdecoderWorkflowConfig.getRequired(task_param_map, CPTMdecoderWorkflowConfig.PARAM_MOD_FILE_PATH, 'modification file path');
            cfg.fixed_mod = CPTMdecoderWorkflowConfig.getRequired(task_param_map, CPTMdecoderWorkflowConfig.PARAM_FIXED_MOD, 'fixed modification list');
            cfg.variable_mod = CPTMdecoderWorkflowConfig.getRequired(task_param_map, CPTMdecoderWorkflowConfig.PARAM_VARIABLE_MOD, 'variable modification list');
            cfg.spec_dir_path = CPTMdecoderWorkflowConfig.getRequired(task_param_map, CPTMdecoderWorkflowConfig.PARAM_SPEC_DIR_PATH, 'spectrum directory path');

            ms1_tol_value = str2double(CPTMdecoderWorkflowConfig.getRequired(task_param_map, CPTMdecoderWorkflowConfig.PARAM_MS1_TOLERANCE_VALUE, 'MS1 tolerance value'));
            ms1_tol_type = strtrim(CPTMdecoderWorkflowConfig.getRequired(task_param_map, CPTMdecoderWorkflowConfig.PARAM_MS1_TOLERANCE_TYPE, 'MS1 tolerance type'));
            cfg.ms1_tolerance = struct('value', ms1_tol_value, 'isppm', strcmpi(ms1_tol_type, 'PPM'));

            cfg.ms2_tolerance = str2double(CPTMdecoderWorkflowConfig.getRequired(task_param_map, CPTMdecoderWorkflowConfig.PARAM_MS2_TOLERANCE, 'MS2 tolerance value'));
            cfg.alpha = str2double(CPTMdecoderWorkflowConfig.getRequired(task_param_map, CPTMdecoderWorkflowConfig.PARAM_ALPHA, 'noise filtering threshold'));
            cfg.fasta_file_path = CPTMdecoderWorkflowConfig.getRequired(task_param_map, CPTMdecoderWorkflowConfig.PARAM_FASTA_FILE_PATH, 'fasta file path');
            cfg.regular_express = CPTMdecoderWorkflowConfig.getRequired(task_param_map, CPTMdecoderWorkflowConfig.PARAM_REGULAR_EXPRESS, 'protein name regex');
            cfg.filtered_res_file_path = CPTMdecoderWorkflowConfig.getOptional(task_param_map, CPTMdecoderWorkflowConfig.PARAM_FILTERED_RES_FILE_PATH, '');
            cfg.model = str2double(CPTMdecoderWorkflowConfig.getRequired(task_param_map, CPTMdecoderWorkflowConfig.PARAM_MODEL, 'algorithm model'));
            cfg.method = str2double(CPTMdecoderWorkflowConfig.getRequired(task_param_map, CPTMdecoderWorkflowConfig.PARAM_METHOD, 'algorithm method'));

            lambda_value = CPTMdecoderWorkflowConfig.getOptional(task_param_map, CPTMdecoderWorkflowConfig.PARAM_LAMBDA, '0');
            cfg.lambda = str2double(lambda_value);
            if isequal(strtrim(CPTMdecoderWorkflowConfig.getOptional(task_param_map, CPTMdecoderWorkflowConfig.PARAM_METHOD, '')), '2') && ...
                    ~task_param_map.isKey(CPTMdecoderWorkflowConfig.PARAM_LAMBDA)
                error('CPTMdecoderWorkflowConfig:MissingLambda', ...
                    'The lasso parameter ''lambda'' is required when method=2.');
            end

            cfg.result_filter_threshold = str2double(CPTMdecoderWorkflowConfig.getRequired(task_param_map, CPTMdecoderWorkflowConfig.PARAM_RESULT_FILTER_THRESHOLD, 'result filter threshold'));
            cfg.enzyme_name = CPTMdecoderWorkflowConfig.getRequired(task_param_map, CPTMdecoderWorkflowConfig.PARAM_ENZYME_NAME, 'enzyme name');
            cfg.enzyme_limits = str2num(CPTMdecoderWorkflowConfig.getRequired(task_param_map, CPTMdecoderWorkflowConfig.PARAM_ENZYME_LIMIT_C_TERM_POSSIBLE_MOD, 'enzyme C-term possible modifications')); %#ok<ST2NM>
            cfg.output_dir_path = CPTMdecoderWorkflowConfig.getRequired(task_param_map, CPTMdecoderWorkflowConfig.PARAM_OUTPUT_DIR_PATH, 'output directory');

            cfg.min_MSMS_num = str2double(CPTMdecoderWorkflowConfig.getOptional(task_param_map, CPTMdecoderWorkflowConfig.PARAM_MIN_MSMS_NUM, '1'));
            cfg.ion_types = [1,2];
            cfg.case_penalty_intens = 'intens_sum';
            cfg.grid_penalty_intens = 'intens_sum';
            cfg.case_OLS_intens_weight = 'none';

            if strcmp(mode, CPTMdecoderWorkflowConfig.ACTION_MSMS_PEPTIDE) || strcmp(mode, CPTMdecoderWorkflowConfig.ACTION_MSMS_ONLY)
                cfg.pep_spec_file_path = CPTMdecoderWorkflowConfig.getRequired(task_param_map, CPTMdecoderWorkflowConfig.PARAM_PEP_SPEC_FILE_PATH, 'peptide-spectrum list path');
                cfg.checked_peptides_res_path = [];
                cfg.msms_res_path = [];
            elseif strcmp(mode, CPTMdecoderWorkflowConfig.ACTION_PEPTIDE_REQUANT)
                cfg.pep_spec_file_path = '';
                cfg.checked_peptides_res_path = CPTMdecoderWorkflowConfig.getRequired(task_param_map, CPTMdecoderWorkflowConfig.PARAM_CHECKED_PEPTIDES_RES_PATH, 'checked peptide result path');
                cfg.msms_res_path = CPTMdecoderWorkflowConfig.getRequired(task_param_map, CPTMdecoderWorkflowConfig.PARAM_MSMS_RES_PATH, 'MSMS result path');
            else
                cfg.pep_spec_file_path = '';
                cfg.checked_peptides_res_path = [];
                cfg.msms_res_path = CPTMdecoderWorkflowConfig.getRequired(task_param_map, CPTMdecoderWorkflowConfig.PARAM_MSMS_RES_PATH, 'MSMS result path');
            end
        end

        function cfg = buildNormalizationMinimalConfigStructFromMap(task_param_map)
            % Build minimal normalization config struct from parameter map.
            %   task_param_map (containers.Map)
            %       parameter key-value map
            % Build minimal normalization config struct from parameter map.
            %   task_param_map (containers.Map)
            %       parameter key-value map
            % Build minimal normalization config consumed by normalization services.
            % This is a compatibility bridge and intentionally keeps only required fields.
            cfg = struct();
            cfg.spec_dir_path = CPTMdecoderWorkflowConfig.getRequired(task_param_map, CPTMdecoderWorkflowConfig.PARAM_SPEC_DIR_PATH, 'spectrum directory path');

            ms1_tol_value = str2double(CPTMdecoderWorkflowConfig.getRequired(task_param_map, CPTMdecoderWorkflowConfig.PARAM_MS1_TOLERANCE_VALUE, 'MS1 tolerance value'));
            ms1_tol_type = strtrim(CPTMdecoderWorkflowConfig.getRequired(task_param_map, CPTMdecoderWorkflowConfig.PARAM_MS1_TOLERANCE_TYPE, 'MS1 tolerance type'));
            cfg.ms1_tolerance = struct('value', ms1_tol_value, 'isppm', strcmpi(ms1_tol_type, 'PPM'));

            cfg.alpha = str2double(CPTMdecoderWorkflowConfig.getRequired(task_param_map, CPTMdecoderWorkflowConfig.PARAM_ALPHA, 'noise filtering threshold'));
            cfg.result_filter_threshold = str2double(CPTMdecoderWorkflowConfig.getRequired(task_param_map, CPTMdecoderWorkflowConfig.PARAM_RESULT_FILTER_THRESHOLD, 'result filter threshold'));
            cfg.output_dir_path = CPTMdecoderWorkflowConfig.getRequired(task_param_map, CPTMdecoderWorkflowConfig.PARAM_OUTPUT_DIR_PATH, 'output directory');
            cfg.min_MSMS_num = str2double(CPTMdecoderWorkflowConfig.getOptional(task_param_map, CPTMdecoderWorkflowConfig.PARAM_MIN_MSMS_NUM, '1'));
            cfg.checked_peptides_res_path = CPTMdecoderWorkflowConfig.getOptional(task_param_map, CPTMdecoderWorkflowConfig.PARAM_CHECKED_PEPTIDES_RES_PATH, []);
        end

        function cfg = buildNormalizationQuantConfigFromMap(task_param_map)
            % Build normalization quant config from parameter map.
            %   task_param_map (containers.Map)
            %       parameter key-value map
            base_cfg_struct = CPTMdecoderWorkflowConfig.buildNormalizationMinimalConfigStructFromMap(task_param_map);
            pair_file_path = CPTMdecoderWorkflowConfig.getRequired(task_param_map, ...
                CPTMdecoderWorkflowConfig.PARAM_NORM_PROTEIN_PEPTIDE_PAIR_PATH, ...
                'normalization protein-peptide pair file path');
            [prot_list, peptide_list] = CProteinPeptidePairFileReader.readPairFile(pair_file_path);

            cfg = struct();
            cfg.msms_cfg = base_cfg_struct;
            cfg.peptide_list = {peptide_list};
            cfg.prot_list = {prot_list};
            cfg.filtered_res_file_path = CPTMdecoderWorkflowConfig.getRequired(task_param_map, CPTMdecoderWorkflowConfig.PARAM_FILTERED_RES_FILE_PATH, 'filtered result path');
            cfg.output_file_name = CPTMdecoderWorkflowConfig.getOptional(task_param_map, CPTMdecoderWorkflowConfig.PARAM_NORM_OUTPUT_FILE_NAME, 'peptide4normalization.txt');
        end

        function cfg = buildNormalizationRequantConfigFromMap(task_param_map)
            % Build normalization requant config from parameter map.
            %   task_param_map (containers.Map)
            %       parameter key-value map
            base_cfg_struct = CPTMdecoderWorkflowConfig.buildNormalizationMinimalConfigStructFromMap(task_param_map);
            cfg = struct();
            cfg.msms_cfg = base_cfg_struct;
            cfg.output_file_name = CPTMdecoderWorkflowConfig.getOptional(task_param_map, CPTMdecoderWorkflowConfig.PARAM_NORM_REQUANT_OUTPUT_FILE_NAME, 'peptide4normalization_requant.txt');
        end

        function cfg = buildSiteConfigStructFromMap(task_param_map)
            % Build site config struct from parameter map.
            %   task_param_map (containers.Map)
            %       parameter key-value map
            cfg = struct();

            output_dir_path = CPTMdecoderWorkflowConfig.getOptional(task_param_map, CPTMdecoderWorkflowConfig.PARAM_OUTPUT_DIR_PATH, '');
            pep_default = fullfile(output_dir_path, 'report_peptide_all.txt');
            intere_default = fullfile(output_dir_path, 'report_site.txt');
            unintere_default = fullfile(output_dir_path, 'report_peptide_uninterested.txt');

            cfg.input_path = CPTMdecoderWorkflowConfig.getOptional(task_param_map, CPTMdecoderWorkflowConfig.PARAM_PEP_LEVEL_FILE_PATH, pep_default);
            if isempty(cfg.input_path)
                cfg.input_path = pep_default;
            end

            cfg.output_intere_path = CPTMdecoderWorkflowConfig.getOptional(task_param_map, CPTMdecoderWorkflowConfig.PARAM_OUTPUT_INTERE_PATH, intere_default);
            if isempty(cfg.output_intere_path)
                cfg.output_intere_path = intere_default;
            end

            cfg.output_unintere_path = CPTMdecoderWorkflowConfig.getOptional(task_param_map, CPTMdecoderWorkflowConfig.PARAM_OUTPUT_UNINTERE_PATH, unintere_default);
            if isempty(cfg.output_unintere_path)
                cfg.output_unintere_path = unintere_default;
            end

            protein_name_abbr_num = str2double(CPTMdecoderWorkflowConfig.getRequired(task_param_map, CPTMdecoderWorkflowConfig.PARAM_PROTEIN_NAME_ABBR_NUM, 'number of protein abbreviation mappings'));
            protein_name_abbr = containers.Map;
            for idx = 1:protein_name_abbr_num
                key_name = [CPTMdecoderWorkflowConfig.PARAM_PREFIX_PROTEIN_NAME_ABBR, num2str(idx)];
                pair_str = CPTMdecoderWorkflowConfig.getRequired(task_param_map, key_name, 'protein abbreviation pair');
                split_str = strsplit(pair_str, '>');
                protein_name_abbr(strtrim(split_str{1})) = strtrim(split_str{2});
            end
            cfg.protein_name_abbr = protein_name_abbr;

            mod_name_abbr_num = str2double(CPTMdecoderWorkflowConfig.getRequired(task_param_map, CPTMdecoderWorkflowConfig.PARAM_MOD_NAME_ABBR_NUM, 'number of modification abbreviation mappings'));
            mod_name_abbr = containers.Map;
            for idx = 1:mod_name_abbr_num
                key_name = [CPTMdecoderWorkflowConfig.PARAM_PREFIX_MOD_NAME_ABBR, num2str(idx)];
                pair_str = CPTMdecoderWorkflowConfig.getRequired(task_param_map, key_name, 'modification abbreviation pair');
                split_str = strsplit(pair_str, '>');
                mod_name_abbr(strtrim(split_str{1})) = strtrim(split_str{2});
            end
            cfg.mod_name_abbr = mod_name_abbr;

            ignore_str = CPTMdecoderWorkflowConfig.getRequired(task_param_map, CPTMdecoderWorkflowConfig.PARAM_IGNORE_STRINGS_SITE_LEVEL, 'ignore strings for site-level summary');
            cfg.ignore_strings = CPTMdecoderWorkflowConfig.parseIgnoreStrings(ignore_str);
            cfg.column_idxs = struct('icol_seq', 2, 'icol_auc', 8);
        end

        function pair_configs = buildMergeToPairConfigsFromMap(task_param_map)
            % Build merge to pair configs from parameter map.
            %   task_param_map (containers.Map)
            %       parameter key-value map
            left_right_out_num = str2double(CPTMdecoderWorkflowConfig.getRequired(task_param_map, CPTMdecoderWorkflowConfig.PARAM_LEFT_RIGHT_OUT_NUM, 'number of pairwise comparisons'));
            left_name = CPTMdecoderWorkflowConfig.getRequired(task_param_map, CPTMdecoderWorkflowConfig.PARAM_LEFT_NAME, 'left group name');
            right_name = CPTMdecoderWorkflowConfig.getRequired(task_param_map, CPTMdecoderWorkflowConfig.PARAM_RIGHT_NAME, 'right group name');
            ignore_str = CPTMdecoderWorkflowConfig.getRequired(task_param_map, CPTMdecoderWorkflowConfig.PARAM_IGNORE_STRINGS_PAIR_LEVEL, 'ignore strings for pair-level merge');
            ignore_strings = CPTMdecoderWorkflowConfig.parseIgnoreStrings(ignore_str);

            pair_configs = cell(left_right_out_num, 1);
            for idx = 1:left_right_out_num
                key_name = [CPTMdecoderWorkflowConfig.PARAM_PREFIX_LEFT_RIGHT_OUT, num2str(idx)];
                left_right_out_str = CPTMdecoderWorkflowConfig.getRequired(task_param_map, key_name, 'pairwise input/output mapping');
                split_str = strsplit(left_right_out_str, {'|', '>'});
                pair_row = strtrim(split_str);
                pair_configs{idx} = CMergeEachPairConfig.fromPairRow(pair_row, left_name, right_name, ignore_strings);
            end
        end

        function cfg = buildMergePairsConfigStructFromMap(task_param_map)
            % Build merge pairs config struct from parameter map.
            %   task_param_map (containers.Map)
            %       parameter key-value map
            pair_num = str2double(CPTMdecoderWorkflowConfig.getRequired(task_param_map, CPTMdecoderWorkflowConfig.PARAM_PAIR_NUM, 'number of pairs to merge'));
            result_paths = cell(pair_num, 1);
            group_titles = cell(pair_num, 2);
            for idx = 1:pair_num
                pair_key = [CPTMdecoderWorkflowConfig.PARAM_PREFIX_PAIR, num2str(idx)];
                lr_key = [CPTMdecoderWorkflowConfig.PARAM_PREFIX_LEFT_RIGHT_NAME, num2str(idx)];
                result_paths{idx} = CPTMdecoderWorkflowConfig.getRequired(task_param_map, pair_key, 'pair-level result path');

                lr_str = CPTMdecoderWorkflowConfig.getRequired(task_param_map, lr_key, 'left-right group names');
                split_lr = strsplit(lr_str, '|');
                group_titles{idx, 1} = strtrim(split_lr{1});
                group_titles{idx, 2} = strtrim(split_lr{2});
            end

            cfg = struct();
            cfg.result_paths = result_paths;
            cfg.output_path = CPTMdecoderWorkflowConfig.getRequired(task_param_map, CPTMdecoderWorkflowConfig.PARAM_FINAL_OUTPUT_PATH, 'final merged output path');
            cfg.group_titles = group_titles;
        end

        function ignore_strings = parseIgnoreStrings(ignore_strings_str)
            % Parse ignore strings from a string.
            %   ignore_strings_str (1 x 1 char/string)
            %       input string containing ignore strings
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
