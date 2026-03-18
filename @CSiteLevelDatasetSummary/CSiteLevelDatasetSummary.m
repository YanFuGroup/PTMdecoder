classdef CSiteLevelDatasetSummary
    % Build site-by-dataset matrix summary from peptide-level result.

    properties
        % Paths
        m_input_path                    % Input path of peptide-level result
        m_output_site_dataset_matrix_path % Output path of site-by-dataset matrix

        % Mapping context
        m_protein_abbr_file_path             % Protein-abbreviation TSV path
        m_protein_abbr_file_col_protein_name % Protein-name column in TSV
        m_protein_abbr_file_col_abbr_name    % Abbreviation column in TSV
        m_protein_name_abbr                  % Protein full name -> abbreviation map
        m_mod_name_abbr                      % Modification full name -> abbreviation map
        m_ignore_strings                     % Strings to remove from modified peptide sequence
        m_column_idxs                        % Column index settings

        % Matrix result
        m_site_names         % 1 x N cell: site names
        m_dataset_names      % 1 x M cell: dataset names
        m_site_dataset_sum   % containers.Map(site -> containers.Map(dataset -> sum))
    end

    methods
        function obj = CSiteLevelDatasetSummary(config)
            % Construct an instance from CSiteLevelDatasetPipelineConfig.
            % Inputs:
            %   config (CSiteLevelDatasetPipelineConfig)
            %       dataset-level site summary configuration
            if nargin < 1 || isempty(config)
                error('CSiteLevelDatasetSummary:MissingConfig', ...
                    'CSiteLevelDatasetPipelineConfig is required.');
            end
            if ~isa(config, 'CSiteLevelDatasetPipelineConfig')
                error('CSiteLevelDatasetSummary:InvalidConfigType', ...
                    'config must be a CSiteLevelDatasetPipelineConfig object.');
            end

            obj.m_input_path = config.input_path;
            obj.m_output_site_dataset_matrix_path = config.output_site_dataset_matrix_path;
            obj.m_protein_abbr_file_path = config.protein_abbr_file_path;
            obj.m_protein_abbr_file_col_protein_name = config.protein_abbr_file_col_protein_name;
            obj.m_protein_abbr_file_col_abbr_name = config.protein_abbr_file_col_abbr_name;
            obj.m_protein_name_abbr = config.protein_name_abbr;
            obj.m_mod_name_abbr = config.mod_name_abbr;
            obj.m_ignore_strings = config.ignore_strings;
            obj.m_column_idxs = config.column_idxs;

            obj.m_site_names = {};
            obj.m_dataset_names = {};
            obj.m_site_dataset_sum = containers.Map('KeyType', 'char', 'ValueType', 'any');
        end


        function run(obj)
            % Build matrix summary and write to output file.
            obj = obj.site_level_dataset_summary();
            obj.write_file();
        end


        % Build matrix summary (scaffold stage).
        obj = site_level_dataset_summary(obj);


        % Write matrix summary to file.
        write_file(obj);
    end
end
