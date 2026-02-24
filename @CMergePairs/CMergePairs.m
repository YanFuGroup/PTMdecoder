classdef CMergePairs
    % Merge pairs and get the joint to form the final result
    
    properties
        % Paths
        m_result_paths;   % The pathes of the paired result
        m_output_path;    % The path of the output

        % Parameters of the merging process
        m_column_idxs;    % The column indexes of the data to be merged
        m_group_titles;   % The titles of the groups

        % The result of the merging
        m_result;         % The result of the merging, which is a cell array
    end
    
    methods
        function obj = CMergePairs(config)
            % Construct an instance of this class
            % Input:
            %   config (CMergePairsConfig)
            %       config object for merge-pairs stage
            %       - result_paths (N x 1 cell)
            %           paths of pair-level merged result files
            %       - output_path (1 x 1 char/string)
            %           output path of final merged report
            %       - group_titles (N x 2 cell)
            %           names of each pair, corresponding to result_paths
            %       - column_idxs (struct)
            %           parser index settings:
            %             * icol_site, icol_pep, icol_charge
            %             * icol_quant_1, icol_quant_2
            if nargin < 1 || isempty(config)
                error('CMergePairs:MissingConfig', ...
                    'CMergePairsConfig is required.');
            end
            if ~isa(config, 'CMergePairsConfig')
                error('CMergePairs:InvalidConfigType', ...
                    'config must be a CMergePairsConfig object.');
            end

            obj.m_result_paths = config.result_paths;
            obj.m_output_path = config.output_path;
            obj.m_group_titles = config.group_titles;
            obj.m_column_idxs = config.column_idxs;
        end
        
        function merge_and_write(obj)
            % Merge the pairs and write the result to the output path

            % Read the results and merge them
            obj.m_result = obj.get_merge_res_of_pairs();

            % Write the result
            obj.write_result();
        end

        % Merge the results of the pairs
        res_joint = get_merge_res_of_pairs(obj);

        % Write the result to the output path
        write_result(obj);
    end
end

