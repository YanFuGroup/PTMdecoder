classdef CMergeEachPairConfig
    % Config for CMergeEachPair.

    properties
        result_path1   % Path of the first site-level result
        result_path2   % Path of the second site-level result
        output_path    % Output path of merged pair report
        group_titles   % 1 x 2 cell: titles of compared groups
        ignore_strings % 1 x M cell: strings ignored when matching peptides
        column_idxs    % struct: column index settings used by parser/merger
    end

    methods
        function obj = CMergeEachPairConfig(cfg)
            % Construct from config struct.
            % Input:
            %   cfg (struct)
            %       configurable fields:
            %       - result_path1, result_path2, output_path
            %       - group_titles, ignore_strings, column_idxs
            if nargin < 1
                cfg = struct();
            end
            cfg = CMergeEachPairConfig.finalize(cfg);

            obj.result_path1 = cfg.result_path1;
            obj.result_path2 = cfg.result_path2;
            obj.output_path = cfg.output_path;
            obj.group_titles = cfg.group_titles;
            obj.ignore_strings = cfg.ignore_strings;
            obj.column_idxs = cfg.column_idxs;
        end
    end

    methods (Static)
        cfgs = fromParamMap(task_param_map)

        function obj = fromPairRow(pair_row, left_name, right_name, ignore_strings)
            % Build one-pair merge config from a [left right output] row.
            % Input:
            %   pair_row (1 x 3 cell)
            %       {left_result_path, right_result_path, output_path}
            %   left_name, right_name (char/string)
            %       names used in report header
            %   ignore_strings (1 x M cell)
            % Output:
            %   obj (CMergeEachPairConfig)
            cfg = struct( ...
                'result_path1', pair_row{1}, ...
                'result_path2', pair_row{2}, ...
                'output_path', pair_row{3}, ...
                'group_titles', {{left_name, right_name}}, ...
                'ignore_strings', {ignore_strings});
            obj = CMergeEachPairConfig(cfg);
        end
    end

    methods (Static, Access = private)
        function cfg = finalize(cfg)
            if ~isfield(cfg, 'result_path1'); cfg.result_path1 = ''; end
            if ~isfield(cfg, 'result_path2'); cfg.result_path2 = ''; end
            if ~isfield(cfg, 'output_path'); cfg.output_path = ''; end
            if ~isfield(cfg, 'group_titles') || isempty(cfg.group_titles); cfg.group_titles = {'left', 'right'}; end
            if ~isfield(cfg, 'ignore_strings') || isempty(cfg.ignore_strings); cfg.ignore_strings = {}; end
            if ~iscell(cfg.ignore_strings); cfg.ignore_strings = {cfg.ignore_strings}; end
            if ~isfield(cfg, 'column_idxs') || isempty(cfg.column_idxs)
                cfg.column_idxs = struct();
                cfg.column_idxs.icol_site = 1;
                cfg.column_idxs.icol_pep = 2;
                cfg.column_idxs.icol_charge = 3;
                cfg.column_idxs.icol_dataset = 4;
                cfg.column_idxs.icol_quant = 8;
                cfg.column_idxs.icol_max = 8;
            end

            if isempty(cfg.result_path1) || isempty(cfg.result_path2) || isempty(cfg.output_path)
                error('CMergeEachPairConfig:InvalidPaths', 'result_path1, result_path2, and output_path must be provided.');
            end
            if ~iscell(cfg.group_titles) || numel(cfg.group_titles) ~= 2
                error('CMergeEachPairConfig:InvalidGroupTitles', 'group_titles must be a 1x2 cell array.');
            end
        end
    end
end
