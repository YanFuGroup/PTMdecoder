classdef CSiteLevelSummary
    % Summary in site level
    
    properties
        % Paths
        m_input_path;               % Input path
        m_output_path_interested;   % Output path for interested sites
        m_output_path_uninterested; % Output path for uninterested sites

        m_protein_abbr_input_mode;               % Protein abbreviation source mode: manual | file

        % Strings
        m_protein_name_abbr;    % Protein name and correspongding abbreviation (map)
        m_mod_name_abbr;        % Modification name and correspongding abbreviation (map)
        m_ignore_strings;       % Ignore strings, using for ignore heavy label modification strings in SILAC data

        % Mapping strategy config
        m_protein_abbr_file_path;                % Protein-abbreviation TSV file path (file mode)
        m_protein_abbr_file_col_protein_name;    % Protein-name column in TSV header (file mode)
        m_protein_abbr_file_col_abbr_name;       % Abbreviation column in TSV header (file mode)

        % Format
        m_column_idxs;          % Column indices of interested sites

        % Site level merging result
        m_result_output_index;          % Site names and corresponding indices of interested sites
        m_result_output_string;         % Output peptide level strings of interested sites
        m_result_output_sum;            % Sum of intensities of interested sites
        m_result_uninterested_string;   % Output peptide level strings of uninterested sites, in peptide level format
    end
    
    methods
        function obj = CSiteLevelSummary(config)
            %CSITELEVELSUMMARY Construct an instance of this class
            % Input:
            %   config (CSiteLevelPipelineConfig)
            %       site-level pipeline configuration
            %       - input_path (1 x 1 char/string)
            %           peptide-level result file path
            %       - output_intere_path (1 x 1 char/string)
            %           output path for interested-site summary
            %       - output_unintere_path (1 x 1 char/string)
            %           output path for uninterested peptide records
            %       - protein_abbr_input_mode (1 x 1 char/string)
            %           protein abbreviation source mode (manual | file)
            %       - protein_abbr_file_path (1 x 1 char/string)
            %           protein-abbreviation TSV path (used in file mode)
            %       - protein_abbr_file_col_protein_name (1 x 1 char/string)
            %           protein-name column in TSV header (used in file mode)
            %       - protein_abbr_file_col_abbr_name (1 x 1 char/string)
            %           abbreviation column in TSV header (used in file mode)
            %       - protein_name_abbr (containers.Map)
            %           mapping: full protein name -> abbreviation
            %       - mod_name_abbr (containers.Map)
            %           mapping: full modification name -> abbreviation
            %       - ignore_strings (1 x M cell)
            %           strings to remove from modified peptide sequence
            %       - column_idxs (struct)
            %           required fields:
            %             * icol_seq : column index of modified sequence
            %             * icol_auc : column index of quantification value
            if nargin < 1 || isempty(config)
                error('CSiteLevelSummary:MissingConfig', ...
                    'CSiteLevelPipelineConfig is required.');
            end
            if ~isa(config, 'CSiteLevelPipelineConfig')
                error('CSiteLevelSummary:InvalidConfigType', ...
                    'config must be a CSiteLevelPipelineConfig object.');
            end

            obj.m_input_path = config.input_path;
            obj.m_output_path_interested = config.output_intere_path;
            obj.m_output_path_uninterested = config.output_unintere_path;
            obj.m_protein_abbr_input_mode = config.protein_abbr_input_mode;
            obj.m_protein_abbr_file_path = config.protein_abbr_file_path;
            obj.m_protein_abbr_file_col_protein_name = config.protein_abbr_file_col_protein_name;
            obj.m_protein_abbr_file_col_abbr_name = config.protein_abbr_file_col_abbr_name;
            obj.m_protein_name_abbr = config.protein_name_abbr;
            obj.m_mod_name_abbr = config.mod_name_abbr;
            obj.m_ignore_strings = config.ignore_strings;
            obj.m_column_idxs = config.column_idxs;
        end
        
        function run(obj)
            % Summary and write the result to files
            obj = obj.site_level_summary();
            obj.write_file();
        end

        % Summary in site level
        obj = site_level_summary(obj);

        % Write to file
        write_file(obj);
    end
end

