classdef CTaskParam
    %CTASKPARAM Class for task-specific parameters

    properties(Access=public)
        % Switches for using various modules
        m_msms_peptide_level_on;
        m_peptide_requant_on;
        m_peptide_only_on;
        m_site_level_on;
        m_merge_to_pair_level_on;
        m_merge_pairs_level_on;

        % Parameter on msms & peptide level processing
        m_mod_file_path;
        m_fixed_mod;
        m_variable_mod;
        m_spec_dir_path;
        m_ms1_tolerance;
        m_ms2_tolerance;
        m_alpha;
        m_fasta_file_path;
        m_regular_express;
        m_pep_spec_file_path;
        m_filtered_res_file_path;
        m_model;
        m_method;
        m_lambda;
        m_result_filter_threshold;
        m_enzyme_name;
        m_enzyme_limits;
        m_output_dir_path;
        m_checked_peptides_res_path;
        m_msms_res_path;
        m_min_MSMS_num;    % Minimum number of MSMS spectra for a peptide to be considered

        % Parameter on site level processing
        m_protein_name_abbr_num;
        m_protein_name_abbr;
        m_mod_name_abbr_num;
        m_mod_name_abbr;
        m_pep_level_file_path;
        m_output_intere_path;
        m_output_unintere_path;
        m_ignore_strings_site_level;

        % Parameter on merge to pair level processing
        m_left_right_out_num;
        m_left_right_out;
        m_left_name;
        m_right_name;
        m_ignore_strings_pair_level;
        
        % Parameter on merge pairs level processing
        m_pair_num;
        m_pair;
        m_left_right_name;
        m_final_output_path;
    end

    methods
        function obj = CTaskParam(param_file_path)
            % Input:
            %   param_file_path (1 x 1 char/string)
            %       task parameter file path

            % Parameter dictionary
            task_param_map = obj.parse_file(param_file_path); 

            % Check if all required parameters are included, and 
            %   throw an error with a specific message if a certain attribute is not found
            obj.check_required_keys(task_param_map);

            % Assign values to the task parameters object
            obj = assign_values(obj, task_param_map);
        end
    end
end

