classdef CSiteLevelSummary
    % Summary in site level
    
    properties
        % Paths
        m_input_path;               % Input path
        m_output_path_interested;   % Output path for interested sites
        m_output_path_uninterested; % Output path for uninterested sites

        % Strings
        m_protein_name_abbr;    % Protein name and correspongding abbreviation (map)
        m_mod_name_abbr;        % Modification name and correspongding abbreviation (map)
        m_ignore_strings;       % Ignore strings, using for ignore heavy label modification strings in SILAC data

        % Format
        m_column_idxs;          % Column indices of interested sites

        % Site level merging result
        m_result_output_index;          % Site names and corresponding indices of interested sites
        m_result_output_string;         % Output peptide level strings of interested sites
        m_result_output_sum;            % Sum of intensities of interested sites
        m_result_uninterested_string;   % Output peptide level strings of uninterested sites, in peptide level format
    end
    
    methods
        function obj = CSiteLevelSummary(input_path_or_task_param_obj, ...
                output_intere_path, output_unintere_path, ...
                protein_name_abbr, mod_name_abbr, ignore_strings, column_idxs)
            %CSITELEVELSUMMARY Construct an instance of this class
            % Input:
            %   input_path_or_task_param_obj (char/string or CTaskParam)
            %       peptide-level result path or task parameter object
            %   output_intere_path (1 x 1 char/string)
            %       output path for interested sites
            %   output_unintere_path (1 x 1 char/string)
            %       output path for uninterested sites
            %   protein_name_abbr (containers.Map)
            %       protein name -> abbreviation
            %   mod_name_abbr (containers.Map)
            %       modification name -> abbreviation
            %   ignore_strings (1 x M cell)
            %       strings to remove from peptide sequence
            %   column_idxs (struct, optional)
            %       column index settings
            if nargin == 1
                % Initialize using task parameter obj
                obj.m_input_path = input_path_or_task_param_obj.m_pep_level_file_path;
                obj.m_output_path_interested = input_path_or_task_param_obj.m_output_intere_path;
                obj.m_output_path_uninterested = input_path_or_task_param_obj.m_output_unintere_path;
                obj.m_protein_name_abbr = input_path_or_task_param_obj.m_protein_name_abbr;
                obj.m_mod_name_abbr = input_path_or_task_param_obj.m_mod_name_abbr;
                obj.m_ignore_strings = input_path_or_task_param_obj.m_ignore_strings_site_level;
                % Column indices
                obj.m_column_idxs.icol_seq = 2;
                obj.m_column_idxs.icol_auc = 8;
            else
                % Initialize using input parameters
                obj.m_input_path = input_path_or_task_param_obj;
                obj.m_output_path_interested = output_intere_path;
                obj.m_output_path_uninterested = output_unintere_path;
                obj.m_protein_name_abbr = protein_name_abbr;
                obj.m_mod_name_abbr = mod_name_abbr;
                obj.m_ignore_strings = ignore_strings;
                obj.m_column_idxs = column_idxs;
            end
        end
        
        function summary_and_write(obj)
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

