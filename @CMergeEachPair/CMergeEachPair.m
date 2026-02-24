classdef CMergeEachPair
    % Merge each pair of the input data and compare the result

    % Attention:
    %   In this experiment, only single raw was analysed in each run, so the
    %   peptides will only appear exclusively in one raw. When there are
    %   fractions (peptides appear in several raw), some filters should be
    %   added in site level. In this step, all of the peptides in different
    %   dataset will be sumed. If one want to compare different dataset, then
    %   several runs which only one dataset is contained should be operated.
    
    properties
        % Paths
        m_result_path1;   % The path of the first result
        m_result_path2;   % The path of the second result
        m_output_path;    % The path of the output

        % Parameters of the merging process
        m_column_idxs;    % The column indexes of the data to be merged
        m_group_titles;   % The titles of the groups
        m_ignore_strings;  % The titles of the columns to be ignored

        % The result of the merging
        m_result;         % The result of the merging
    end
    
    methods
        function obj = CMergeEachPair(config)
            % Construct an instance of this class
            % Input:
            %   config (CMergeEachPairConfig)
            %       config object for one pair merge
            %       - result_path1 (1 x 1 char/string)
            %           first site-level result file path
            %       - result_path2 (1 x 1 char/string)
            %           second site-level result file path
            %       - output_path (1 x 1 char/string)
            %           output path for merged pair result
            %       - group_titles (1 x 2 cell)
            %           names used in output header, e.g. {left_name, right_name}
            %       - ignore_strings (1 x M cell)
            %           strings removed from peptide sequence before matching
            %       - column_idxs (struct)
            %           fields used by parser and merger:
            %             * icol_site, icol_pep, icol_charge
            %             * icol_dataset, icol_quant, icol_max
            if nargin < 1 || isempty(config)
                error('CMergeEachPair:MissingConfig', ...
                    'CMergeEachPairConfig is required.');
            end
            if ~isa(config, 'CMergeEachPairConfig')
                error('CMergeEachPair:InvalidConfigType', ...
                    'config must be a CMergeEachPairConfig object.');
            end

            obj.m_result_path1 = config.result_path1;
            obj.m_result_path2 = config.result_path2;
            obj.m_output_path = config.output_path;
            obj.m_group_titles = config.group_titles;
            obj.m_ignore_strings = config.ignore_strings;
            obj.m_column_idxs = config.column_idxs;
            obj.m_result = [];
        end
        
        function merge_and_write(obj)
            % Merge the two results and write the result

            % Read the two results and merge them
            obj.m_result = obj.get_merge_res_of_one_pair();

            % Write the result
            obj.write_result();
        end

        % Get the result of merging one pair of the input data
        res_intersec = get_merge_res_of_one_pair(obj);
        
        % Write the result
        write_result(obj);
    end
end

