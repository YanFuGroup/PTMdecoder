classdef CSiteLevelPipelineConfig
    % Config for CSiteLevelSummary.

    properties
        input_path           % Path of peptide-level result file
        output_intere_path   % Output path for interested-site summary
        output_unintere_path % Output path for uninterested peptide records
        protein_name_abbr    % containers.Map: protein full name -> abbreviation
        mod_name_abbr        % containers.Map: modification full name -> abbreviation
        ignore_strings       % 1 x M cell: strings ignored in modified peptide sequence
        column_idxs          % struct: column index settings (icol_seq, icol_auc)
    end

    methods
        function obj = CSiteLevelPipelineConfig(cfg)
            % Construct from config struct.
            % Input:
            %   cfg (struct)
            %       configurable fields:
            %       - input_path, output_intere_path, output_unintere_path
            %       - protein_name_abbr, mod_name_abbr
            %       - ignore_strings, column_idxs
            if nargin < 1
                cfg = struct();
            end
            cfg = CSiteLevelPipelineConfig.finalize(cfg);

            obj.input_path = cfg.input_path;
            obj.output_intere_path = cfg.output_intere_path;
            obj.output_unintere_path = cfg.output_unintere_path;
            obj.protein_name_abbr = cfg.protein_name_abbr;
            obj.mod_name_abbr = cfg.mod_name_abbr;
            obj.ignore_strings = cfg.ignore_strings;
            obj.column_idxs = cfg.column_idxs;
        end
    end

    methods (Static)
        cfg = fromParamMap(task_param_map)
    end

    methods (Static, Access = private)
        function cfg = finalize(cfg)
            if ~isfield(cfg, 'input_path'); cfg.input_path = ''; end
            if ~isfield(cfg, 'output_intere_path'); cfg.output_intere_path = ''; end
            if ~isfield(cfg, 'output_unintere_path'); cfg.output_unintere_path = ''; end
            if ~isfield(cfg, 'protein_name_abbr') || isempty(cfg.protein_name_abbr); cfg.protein_name_abbr = containers.Map; end
            if ~isfield(cfg, 'mod_name_abbr') || isempty(cfg.mod_name_abbr); cfg.mod_name_abbr = containers.Map; end
            if ~isfield(cfg, 'ignore_strings') || isempty(cfg.ignore_strings); cfg.ignore_strings = {}; end
            if ~iscell(cfg.ignore_strings); cfg.ignore_strings = {cfg.ignore_strings}; end
            if ~isfield(cfg, 'column_idxs') || isempty(cfg.column_idxs)
                cfg.column_idxs = struct('icol_seq', 2, 'icol_auc', 8);
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
