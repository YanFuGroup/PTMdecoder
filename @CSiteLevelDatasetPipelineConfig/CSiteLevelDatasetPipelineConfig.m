classdef CSiteLevelDatasetPipelineConfig
    % Config for site-level dataset matrix summary.

    properties
        input_path                           % Path of peptide-level result file
        output_site_dataset_matrix_path      % Output path of site x dataset matrix
        protein_abbr_file_path               % Path of protein-abbreviation TSV file
        protein_abbr_file_col_protein_name   % Column name for protein name in TSV header
        protein_abbr_file_col_abbr_name      % Column name for abbreviation in TSV header
        protein_name_extract_regex           % Regex pattern for extracting protein key before map lookup
        protein_name_abbr                    % containers.Map: protein full name -> abbreviation
        mod_name_abbr                        % containers.Map: modification full name -> abbreviation
        ignore_strings                       % 1 x M cell: strings ignored in modified peptide sequence
        site_position_count_initial_m        % logical: whether the first initial M is counted in site numbering
        column_idxs                          % struct: column index settings (icol_seq, icol_dataset, icol_auc)
    end

    methods
        function obj = CSiteLevelDatasetPipelineConfig(cfg)
            % Build config object for site-level dataset matrix summary.
            % Inputs:
            %   cfg (struct)
            %       Configurable fields:
            %       - input_path, output_site_dataset_matrix_path
            %       - protein_abbr_file_path, protein_abbr_file_col_protein_name
            %       - protein_abbr_file_col_abbr_name
            %       - protein_name_extract_regex
            %       - protein_name_abbr, mod_name_abbr
            %       - ignore_strings, column_idxs
            if nargin < 1
                cfg = struct();
            end

            cfg = CSiteLevelDatasetPipelineConfig.finalize(cfg);

            obj.input_path = cfg.input_path;
            obj.output_site_dataset_matrix_path = cfg.output_site_dataset_matrix_path;
            obj.protein_abbr_file_path = cfg.protein_abbr_file_path;
            obj.protein_abbr_file_col_protein_name = cfg.protein_abbr_file_col_protein_name;
            obj.protein_abbr_file_col_abbr_name = cfg.protein_abbr_file_col_abbr_name;
            obj.protein_name_extract_regex = cfg.protein_name_extract_regex;
            obj.protein_name_abbr = cfg.protein_name_abbr;
            obj.mod_name_abbr = cfg.mod_name_abbr;
            obj.ignore_strings = cfg.ignore_strings;
            obj.site_position_count_initial_m = cfg.site_position_count_initial_m;
            obj.column_idxs = cfg.column_idxs;
        end
    end

    methods (Static)
        cfg = fromParamMap(task_param_map)
    end

    methods (Static, Access = private)
        function cfg = finalize(cfg)
            if ~isfield(cfg, 'input_path'); cfg.input_path = ''; end
            if ~isfield(cfg, 'output_site_dataset_matrix_path'); cfg.output_site_dataset_matrix_path = ''; end
            if ~isfield(cfg, 'protein_abbr_file_path'); cfg.protein_abbr_file_path = ''; end
            if ~isfield(cfg, 'protein_abbr_file_col_protein_name'); cfg.protein_abbr_file_col_protein_name = ''; end
            if ~isfield(cfg, 'protein_abbr_file_col_abbr_name'); cfg.protein_abbr_file_col_abbr_name = ''; end
            if ~isfield(cfg, 'protein_name_extract_regex'); cfg.protein_name_extract_regex = ''; end
            if ~isfield(cfg, 'protein_name_abbr') || isempty(cfg.protein_name_abbr); cfg.protein_name_abbr = containers.Map; end
            if ~isfield(cfg, 'mod_name_abbr') || isempty(cfg.mod_name_abbr); cfg.mod_name_abbr = containers.Map; end
            if ~isfield(cfg, 'ignore_strings') || isempty(cfg.ignore_strings); cfg.ignore_strings = {}; end
            if ~iscell(cfg.ignore_strings); cfg.ignore_strings = {cfg.ignore_strings}; end
            if ~isfield(cfg, 'site_position_count_initial_m') || isempty(cfg.site_position_count_initial_m)
                cfg.site_position_count_initial_m = true;
            end
            if isstring(cfg.protein_name_extract_regex)
                cfg.protein_name_extract_regex = char(cfg.protein_name_extract_regex);
            end
            cfg.protein_name_extract_regex = strtrim(cfg.protein_name_extract_regex);
            if ~isfield(cfg, 'column_idxs') || isempty(cfg.column_idxs)
                cfg.column_idxs = struct('icol_seq', 2, 'icol_dataset', 4, 'icol_auc', 8);
            end

            if ~isa(cfg.protein_name_abbr, 'containers.Map')
                error('CSiteLevelDatasetPipelineConfig:InvalidProteinNameAbbrMap', ...
                    'protein_name_abbr must be a containers.Map.');
            end
            if ~isa(cfg.mod_name_abbr, 'containers.Map')
                error('CSiteLevelDatasetPipelineConfig:InvalidModNameAbbrMap', ...
                    'mod_name_abbr must be a containers.Map.');
            end

            if isempty(cfg.input_path)
                error('CSiteLevelDatasetPipelineConfig:InvalidInputPath', ...
                    'input_path must be provided.');
            end
            if isempty(cfg.output_site_dataset_matrix_path)
                error('CSiteLevelDatasetPipelineConfig:InvalidOutputPath', ...
                    'output_site_dataset_matrix_path must be provided.');
            end
            if isempty(cfg.protein_abbr_file_path)
                error('CSiteLevelDatasetPipelineConfig:InvalidProteinAbbrFilePath', ...
                    'protein_abbr_file_path must be provided.');
            end
            if isempty(cfg.protein_abbr_file_col_protein_name)
                error('CSiteLevelDatasetPipelineConfig:InvalidProteinNameColumn', ...
                    'protein_abbr_file_col_protein_name must be provided.');
            end
            if isempty(cfg.protein_abbr_file_col_abbr_name)
                error('CSiteLevelDatasetPipelineConfig:InvalidAbbrColumn', ...
                    'protein_abbr_file_col_abbr_name must be provided.');
            end

            if ~isfield(cfg.column_idxs, 'icol_seq') || ...
                    ~isfield(cfg.column_idxs, 'icol_dataset') || ...
                    ~isfield(cfg.column_idxs, 'icol_auc')
                error('CSiteLevelDatasetPipelineConfig:InvalidColumnIdxs', ...
                    'column_idxs must include icol_seq, icol_dataset, and icol_auc.');
            end
        end
    end
end
