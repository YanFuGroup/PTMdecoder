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
        function obj = CMergePairs(result_paths_or_task_param_obj,...
                output_path, group_titles, column_idxs)
            % Construct an instance of this class
            % Input:
            %   result_paths_or_task_param_obj (cell or CTaskParam)
            %       pair result paths or task parameter object
            %   output_path (1 x 1 char/string)
            %       output file path
            %   group_titles (N x 2 cell)
            %       group names for output header
            %   column_idxs (struct, optional)
            %       indices of columns in input files
            if nargin == 1
                task_param = result_paths_or_task_param_obj;
                obj.m_result_paths = task_param.m_pair;
                obj.m_output_path = task_param.m_final_output_path;
                obj.m_group_titles = task_param.m_left_right_name;
            else
                obj.m_result_paths = result_paths_or_task_param_obj;
                obj.m_output_path = output_path;
                obj.m_group_titles = group_titles;
            end
            if nargin < 4
                % Format of lines
                column_idxs.icol_site = 1;
                column_idxs.icol_pep = 2;
                column_idxs.icol_charge = 3;
                column_idxs.icol_quant_1 = 4;
                column_idxs.icol_quant_2 = 5;
            end
            obj.m_column_idxs = column_idxs;
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

