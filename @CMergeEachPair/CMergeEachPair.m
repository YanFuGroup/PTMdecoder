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
        function obj = CMergeEachPair(result_path1,result_path2,output_path,...
                group_titles,ignore_strings,column_idxs)
            % Construct an instance of this class
            obj.m_result_path1 = result_path1;
            obj.m_result_path2 = result_path2;
            obj.m_output_path = output_path;
            obj.m_group_titles = group_titles;
            obj.m_ignore_strings = ignore_strings;
            if nargin < 6
                column_idxs.icol_site = 1;   % The index of column of site name in site lines
                column_idxs.icol_pep = 2;    % The index of column of peptide sequence in peptide lines
                column_idxs.icol_charge = 3; % The index of column of charge in peptide lines
                column_idxs.icol_dataset = 4;% The index of column of dataset name in peptide lines
                column_idxs.icol_quant = 8;  % The index of column of quantification in peptide lines
                column_idxs.icol_max = 8;    % The maximum index of column
            end
            obj.m_column_idxs = column_idxs;
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

