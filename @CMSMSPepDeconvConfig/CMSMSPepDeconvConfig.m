classdef CMSMSPepDeconvConfig
    % Config for CMSMSPepDeconv.

    properties
        mod_file_path
        fixed_mod
        variable_mod
        spec_dir_path
        ms1_tolerance
        ms2_tolerance
        alpha
        fasta_file_path
        regular_express
        pep_spec_file_path
        filtered_res_file_path
        model
        method
        lambda
        result_filter_threshold
        enzyme_name
        enzyme_limits
        output_dir_path
        checked_peptides_res_path
        msms_res_path
        min_MSMS_num
        ion_types
        case_penalty_intens
        grid_penalty_intens
        case_OLS_intens_weight
    end

    methods
        function obj = CMSMSPepDeconvConfig(cfg)
            if nargin < 1
                cfg = struct();
            end
            cfg = CMSMSPepDeconvConfig.finalize(cfg);

            obj.mod_file_path = cfg.mod_file_path;
            obj.fixed_mod = cfg.fixed_mod;
            obj.variable_mod = cfg.variable_mod;
            obj.spec_dir_path = cfg.spec_dir_path;
            obj.ms1_tolerance = cfg.ms1_tolerance;
            obj.ms2_tolerance = cfg.ms2_tolerance;
            obj.alpha = cfg.alpha;
            obj.fasta_file_path = cfg.fasta_file_path;
            obj.regular_express = cfg.regular_express;
            obj.pep_spec_file_path = cfg.pep_spec_file_path;
            obj.filtered_res_file_path = cfg.filtered_res_file_path;
            obj.model = cfg.model;
            obj.method = cfg.method;
            obj.lambda = cfg.lambda;
            obj.result_filter_threshold = cfg.result_filter_threshold;
            obj.enzyme_name = cfg.enzyme_name;
            obj.enzyme_limits = cfg.enzyme_limits;
            obj.output_dir_path = cfg.output_dir_path;
            obj.checked_peptides_res_path = cfg.checked_peptides_res_path;
            obj.msms_res_path = cfg.msms_res_path;
            obj.min_MSMS_num = cfg.min_MSMS_num;
            obj.ion_types = cfg.ion_types;
            obj.case_penalty_intens = cfg.case_penalty_intens;
            obj.grid_penalty_intens = cfg.grid_penalty_intens;
            obj.case_OLS_intens_weight = cfg.case_OLS_intens_weight;
        end
    end

    methods (Static, Access = private)
        function cfg = finalize(cfg)
            if ~isfield(cfg, 'mod_file_path'); cfg.mod_file_path = ''; end
            if ~isfield(cfg, 'fixed_mod'); cfg.fixed_mod = ''; end
            if ~isfield(cfg, 'variable_mod'); cfg.variable_mod = ''; end
            if ~isfield(cfg, 'spec_dir_path'); cfg.spec_dir_path = ''; end
            if ~isfield(cfg, 'ms1_tolerance'); cfg.ms1_tolerance = []; end
            if ~isfield(cfg, 'ms2_tolerance'); cfg.ms2_tolerance = []; end
            if ~isfield(cfg, 'alpha'); cfg.alpha = 0; end
            if ~isfield(cfg, 'fasta_file_path'); cfg.fasta_file_path = ''; end
            if ~isfield(cfg, 'regular_express'); cfg.regular_express = ''; end
            if ~isfield(cfg, 'pep_spec_file_path'); cfg.pep_spec_file_path = ''; end
            if ~isfield(cfg, 'filtered_res_file_path'); cfg.filtered_res_file_path = ''; end
            if ~isfield(cfg, 'model'); cfg.model = []; end
            if ~isfield(cfg, 'method'); cfg.method = []; end
            if ~isfield(cfg, 'lambda') || isempty(cfg.lambda); cfg.lambda = 0; end
            if ~isfield(cfg, 'result_filter_threshold') || isempty(cfg.result_filter_threshold); cfg.result_filter_threshold = 0; end
            if ~isfield(cfg, 'enzyme_name'); cfg.enzyme_name = ''; end
            if ~isfield(cfg, 'enzyme_limits') || isempty(cfg.enzyme_limits); cfg.enzyme_limits = []; end
            if ~isfield(cfg, 'output_dir_path'); cfg.output_dir_path = ''; end
            if ~isfield(cfg, 'checked_peptides_res_path'); cfg.checked_peptides_res_path = []; end
            if ~isfield(cfg, 'msms_res_path'); cfg.msms_res_path = []; end
            if ~isfield(cfg, 'min_MSMS_num') || isempty(cfg.min_MSMS_num); cfg.min_MSMS_num = 1; end
            if ~isfield(cfg, 'ion_types') || isempty(cfg.ion_types); cfg.ion_types = [1,2]; end
            if ~isfield(cfg, 'case_penalty_intens') || isempty(cfg.case_penalty_intens); cfg.case_penalty_intens = 'intens_sum'; end
            if ~isfield(cfg, 'grid_penalty_intens') || isempty(cfg.grid_penalty_intens); cfg.grid_penalty_intens = 'intens_sum'; end
            if ~isfield(cfg, 'case_OLS_intens_weight') || isempty(cfg.case_OLS_intens_weight); cfg.case_OLS_intens_weight = 'none'; end

            if isempty(cfg.mod_file_path)
                error('CMSMSPepDeconvConfig:InvalidModFilePath', 'mod_file_path must be provided.');
            end
            if isempty(cfg.spec_dir_path)
                error('CMSMSPepDeconvConfig:InvalidSpecDirPath', 'spec_dir_path must be provided.');
            end
            if isempty(cfg.fasta_file_path)
                error('CMSMSPepDeconvConfig:InvalidFastaPath', 'fasta_file_path must be provided.');
            end
            if isempty(cfg.output_dir_path)
                error('CMSMSPepDeconvConfig:InvalidOutputDir', 'output_dir_path must be provided.');
            end
        end
    end
end
