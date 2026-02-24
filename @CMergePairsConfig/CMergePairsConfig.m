classdef CMergePairsConfig
    % Config for CMergePairs.

    properties
        result_paths % N x 1 cell: pair-level result file paths
        output_path  % Output path of final merged report
        group_titles % N x 2 cell: pair names corresponding to result_paths
        column_idxs  % struct: column index settings used by parser
    end

    methods
        function obj = CMergePairsConfig(cfg)
            % Construct from config struct.
            % Input:
            %   cfg (struct)
            %       configurable fields:
            %       - result_paths, output_path, group_titles, column_idxs
            if nargin < 1
                cfg = struct();
            end
            cfg = CMergePairsConfig.finalize(cfg);

            obj.result_paths = cfg.result_paths;
            obj.output_path = cfg.output_path;
            obj.group_titles = cfg.group_titles;
            obj.column_idxs = cfg.column_idxs;
        end
    end

    methods (Static, Access = private)
        function cfg = finalize(cfg)
            if ~isfield(cfg, 'result_paths') || isempty(cfg.result_paths); cfg.result_paths = {}; end
            if ~isfield(cfg, 'output_path'); cfg.output_path = ''; end
            if ~isfield(cfg, 'group_titles') || isempty(cfg.group_titles); cfg.group_titles = {}; end
            if ~isfield(cfg, 'column_idxs') || isempty(cfg.column_idxs)
                cfg.column_idxs = struct();
                cfg.column_idxs.icol_site = 1;
                cfg.column_idxs.icol_pep = 2;
                cfg.column_idxs.icol_charge = 3;
                cfg.column_idxs.icol_quant_1 = 4;
                cfg.column_idxs.icol_quant_2 = 5;
            end

            if ~iscell(cfg.result_paths)
                error('CMergePairsConfig:InvalidResultPaths', 'result_paths must be a cell array.');
            end
            if isempty(cfg.output_path)
                error('CMergePairsConfig:InvalidOutputPath', 'output_path must be provided.');
            end
            if ~iscell(cfg.group_titles)
                error('CMergePairsConfig:InvalidGroupTitles', 'group_titles must be a cell array.');
            end
        end
    end
end
