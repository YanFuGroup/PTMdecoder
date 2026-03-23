classdef CSiteLevelPipelineConfig
    % Config for CSiteLevelSummary.

    properties
        input_path           % Path of peptide-level result file
        output_intere_path   % Output path for interested-site summary
        output_unintere_path % Output path for uninterested peptide records
        protein_abbr_input_mode % char: protein abbreviation source mode (manual | file)
        protein_abbr_file_path % Path of protein-abbreviation TSV file in file mode
        protein_abbr_file_col_protein_name % Column name for protein name in TSV header
        protein_abbr_file_col_abbr_name % Column name for abbreviation in TSV header
        protein_name_abbr    % containers.Map: protein full name -> abbreviation
        mod_name_abbr        % containers.Map: modification full name -> abbreviation
        ignore_strings       % 1 x M cell: strings ignored in modified peptide sequence
        site_position_count_initial_m % logical: whether the first initial M is counted in site numbering
        column_idxs          % struct: column index settings (icol_seq, icol_auc)
    end

    methods
        function obj = CSiteLevelPipelineConfig(cfg)
            % Construct from config struct.
            % Input:
            %   cfg (struct)
            %       configurable fields:
            %       - input_path, output_intere_path, output_unintere_path
            %       - protein_abbr_input_mode
            %       - protein_abbr_file_path, protein_abbr_file_col_protein_name
            %       - protein_abbr_file_col_abbr_name
            %       - protein_name_abbr, mod_name_abbr
            %       - ignore_strings, column_idxs
            if nargin < 1
                cfg = struct();
            end
            cfg = CSiteLevelPipelineConfig.finalize(cfg);

            obj.input_path = cfg.input_path;
            obj.output_intere_path = cfg.output_intere_path;
            obj.output_unintere_path = cfg.output_unintere_path;
            obj.protein_abbr_input_mode = cfg.protein_abbr_input_mode;
            obj.protein_abbr_file_path = cfg.protein_abbr_file_path;
            obj.protein_abbr_file_col_protein_name = cfg.protein_abbr_file_col_protein_name;
            obj.protein_abbr_file_col_abbr_name = cfg.protein_abbr_file_col_abbr_name;
            obj.protein_name_abbr = cfg.protein_name_abbr;
            obj.mod_name_abbr = cfg.mod_name_abbr;
            obj.ignore_strings = cfg.ignore_strings;
            obj.site_position_count_initial_m = cfg.site_position_count_initial_m;
            obj.column_idxs = cfg.column_idxs;
        end
    end

    methods (Static)
        cfg = fromParamMap(task_param_map)
        protein_name_abbr_map = loadProteinAbbrMapFromTsv(file_path, col_protein_name, col_abbr_name)
    end

    methods (Static, Access = private)
        function cfg = finalize(cfg)
            if ~isfield(cfg, 'input_path'); cfg.input_path = ''; end
            if ~isfield(cfg, 'output_intere_path'); cfg.output_intere_path = ''; end
            if ~isfield(cfg, 'output_unintere_path'); cfg.output_unintere_path = ''; end
            if ~isfield(cfg, 'protein_abbr_input_mode') || isempty(cfg.protein_abbr_input_mode); cfg.protein_abbr_input_mode = 'manual'; end
            if ~isfield(cfg, 'protein_abbr_file_path'); cfg.protein_abbr_file_path = ''; end
            if ~isfield(cfg, 'protein_abbr_file_col_protein_name'); cfg.protein_abbr_file_col_protein_name = ''; end
            if ~isfield(cfg, 'protein_abbr_file_col_abbr_name'); cfg.protein_abbr_file_col_abbr_name = ''; end
            if ~isfield(cfg, 'protein_name_abbr') || isempty(cfg.protein_name_abbr); cfg.protein_name_abbr = containers.Map; end
            if ~isfield(cfg, 'mod_name_abbr') || isempty(cfg.mod_name_abbr); cfg.mod_name_abbr = containers.Map; end
            if ~isfield(cfg, 'ignore_strings') || isempty(cfg.ignore_strings); cfg.ignore_strings = {}; end
            if ~iscell(cfg.ignore_strings); cfg.ignore_strings = {cfg.ignore_strings}; end
            if ~isfield(cfg, 'site_position_count_initial_m') || isempty(cfg.site_position_count_initial_m)
                cfg.site_position_count_initial_m = false;
            end
            if ~isfield(cfg, 'column_idxs') || isempty(cfg.column_idxs)
                cfg.column_idxs = struct('icol_seq', 2, 'icol_auc', 8);
            end

            if isstring(cfg.protein_abbr_input_mode)
                cfg.protein_abbr_input_mode = char(cfg.protein_abbr_input_mode);
            end
            cfg.protein_abbr_input_mode = lower(strtrim(cfg.protein_abbr_input_mode));

            if ~ismember(cfg.protein_abbr_input_mode, {'manual', 'file'})
                error('CSiteLevelPipelineConfig:InvalidProteinAbbrInputMode', ...
                    'protein_abbr_input_mode must be either ''manual'' or ''file''.');
            end

            if strcmp(cfg.protein_abbr_input_mode, 'file')
                if isempty(cfg.protein_abbr_file_path)
                    error('CSiteLevelPipelineConfig:InvalidProteinAbbrFilePath', ...
                        'protein_abbr_file_path must be provided when protein_abbr_input_mode is ''file''.');
                end
                if isempty(cfg.protein_abbr_file_col_protein_name)
                    error('CSiteLevelPipelineConfig:InvalidProteinAbbrFileProteinNameColumn', ...
                        'protein_abbr_file_col_protein_name must be provided when protein_abbr_input_mode is ''file''.');
                end
                if isempty(cfg.protein_abbr_file_col_abbr_name)
                    error('CSiteLevelPipelineConfig:InvalidProteinAbbrFileAbbrColumn', ...
                        'protein_abbr_file_col_abbr_name must be provided when protein_abbr_input_mode is ''file''.');
                end
            end

            if isempty(cfg.input_path)
                error('CSiteLevelPipelineConfig:InvalidInputPath', 'input_path must be provided.');
            end
            if isempty(cfg.output_intere_path)
                error('CSiteLevelPipelineConfig:InvalidInterestedOutputPath', 'output_intere_path must be provided.');
            end
            if isempty(cfg.output_unintere_path)
                error('CSiteLevelPipelineConfig:InvalidUninterestedOutputPath', 'output_unintere_path must be provided.');
            end
            if ~isfield(cfg.column_idxs, 'icol_seq') || ~isfield(cfg.column_idxs, 'icol_auc')
                error('CSiteLevelPipelineConfig:InvalidColumnIdxs', 'column_idxs must include icol_seq and icol_auc.');
            end
        end
    end
end
