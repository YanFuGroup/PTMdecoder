function main_processes(~, varargin)
% Load the files and run PTMDecoder iteratively
% Input:
%   varargin (1 x N cell)
%       list of task parameter file paths

for i = 1:length(varargin)
    fprintf('Processing %s\n', varargin{i});
    taskParam = CTaskParam(varargin{i});
    
    % Run the msms-peptide-level process
    if taskParam.m_msms_peptide_level_on
        msms_pep_level_process = CMSMSPepDeconv(taskParam);
        msms_pep_level_process.startRun();

    elseif taskParam.m_peptide_requant_on
        requant_process = CMSMSPepDeconv(taskParam);
        requant_process.requant();

    elseif taskParam.m_peptide_only_on
        peptide_only_process = CMSMSPepDeconv(taskParam);
        peptide_only_process.PepLevelRun();
    end

    % The following processes can be run independently
    
    % Run the site-level process
    if taskParam.m_site_level_on
        site_level_process = CSiteLevelSummary(taskParam);
        site_level_process.summary_and_write();
    end

    % Run the merge-to-pair-level process
    if taskParam.m_merge_to_pair_level_on
        for idx_pairs = 1:taskParam.m_left_right_out_num
            current_pair = taskParam.m_left_right_out(idx_pairs, :);
            merge_to_pair_level_process = CMergeEachPair(current_pair{1},...
                current_pair{2}, current_pair{3}, {taskParam.m_left_name, taskParam.m_right_name},...
                taskParam.m_ignore_strings_pair_level);
            merge_to_pair_level_process.merge_and_write();
        end
    end

    % Run the merge-pairs-level process
    if taskParam.m_merge_pairs_level_on
        merge_pairs_level_process = CMergePairs(taskParam);
        merge_pairs_level_process.merge_and_write();
    end
end
end