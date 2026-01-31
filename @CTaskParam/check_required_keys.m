function check_required_keys(~, task_param_map)
% Check if all the required keys are present in the task_param_map.
% Input:
%   task_param_map (containers.Map)
%       task parameter map

% Note: This function assumes that the task_param_map is a containers.Map object.

%% Overall check

% Check if the task_param_map is a containers.Map object.
if ~isa(task_param_map, 'containers.Map')
    error('The task_param_map should be a containers.Map object.')
end

% Check if the task_param_map contains the keys of which modules are used
% If the key is present and its value is '1', set the corresponding obj.m_ variable to true, 
%   otherwise set it to false.
module_names = { ...
    'msms_peptide_level_on', ...
    'peptide_requant_on', ...
    'peptide_only_on', ...          % This is for the top-1 result of search engine
    'site_level_on', ...
    'merge_to_pair_level_on', ...
    'merge_pairs_level_on' ...
    };

% Check the runs of peptides, can't at the same time
count = 0;
for i = 1:3
    if task_param_map.isKey(module_names{i}) && strcmp(task_param_map(module_names{i}), '1')
        count = count + 1;
    end
end
if count > 1
    error('Cannot process the peptide level quantifications at the same time.')
end

%% MSMS & Peptide level processing
if task_param_map.isKey(module_names{1}) && strcmp(task_param_map(module_names{1}), '1')
    check_required_keys_msms_peptide_level(task_param_map)

%% Peptide requant processing
elseif task_param_map.isKey(module_names{2}) && strcmp(task_param_map(module_names{2}), '1')
    check_required_keys_peptide_requant(task_param_map)

%% Peptide only processing
elseif task_param_map.isKey(module_names{3}) && strcmp(task_param_map(module_names{3}), '1')
    check_required_keys_peptide_only(task_param_map)

end

%% Site level processing
if task_param_map.isKey(module_names{4}) && strcmp(task_param_map(module_names{4}), '1')
    check_required_keys_site_level(task_param_map)
end

%% Merge to pair level processing
if task_param_map.isKey(module_names{5}) && strcmp(task_param_map(module_names{5}), '1')
    check_required_keys_merge_to_pair_level(task_param_map)
end

%% Merge pairs level processing
if task_param_map.isKey(module_names{6}) && strcmp(task_param_map(module_names{6}), '1')
    check_required_keys_merge_pairs_level(task_param_map)
end
end



function check_required_keys_msms_peptide_level(task_param_map)
% Check if all the required keys are present in the task_param_map for MSMS & Peptide level processing.
% Input:
%   task_param_map (containers.Map)
%       task parameter map
    
key_needed_name_explain = { ...
    % Parameter on msms & peptide level processing
    'mod_file_path', 'the path of the modification database file';...
    'fixed_mod', 'the fixed modifications';...
    'variable_mod', 'the variable modifications';...
    'spec_dir_path', 'the path of directory containing the mgf file(s)';...
    'ms1_tolerance_value', 'the value of MS1 tolerance';...
    'ms1_tolerance_type', 'the type of MS1 tolerance';...
    'ms2_tolerance', 'the value of MS2 tolerance';...
    'alpha', 'the relative intensity threshold using in denoising';...
    'fasta_file_path', 'the path of the fasta file';...
    'regular_express', 'the regular expression for parsing the protein names';...
    'pep_spec_file_path', 'the path of the peptide-spectra list file';...
    'model', 'the model used in algorithm (FEV/FEC/FEE)';...
    'method', 'the method used in algorithm (OLS/lasso)';...
    'result_filter_threshold', 'the result filter threshold of relative abundance';...
    'enzyme_name', 'the name of enzyme';...
    'enzyme_limit_C_term_possible_mod', 'the feasible C-terminal modification defined by the enzyme';...
    'output_dir_path', 'the path of output directory'
};

% Check if all the required keys are present in the task_param_map
check_keys_in_map(task_param_map, key_needed_name_explain);

% In core model, lambda is used in lasso model.
% If the 'method' key is set to '2', the function checks if the 'lambda' key is present.
if isequal(task_param_map('method'),'2') && ~task_param_map.isKey('lambda')
    error('The param lasso parameter ''lambda'' is not found.')
end
end



function check_required_keys_peptide_requant(task_param_map)
% Check if all the required keys are present in the task_param_map for MSMS & Peptide level processing.
% Input:
%   task_param_map (containers.Map)
%       task parameter map

% TODO: find the minimal set of required keys for peptide level requantification

key_needed_name_explain = { ...
    % Parameter on msms & peptide level processing
    'mod_file_path', 'the path of the modification database file';...
    'fixed_mod', 'the fixed modifications';...
    'variable_mod', 'the variable modifications';...
    'spec_dir_path', 'the path of directory containing the mgf file(s)';...
    'ms1_tolerance_value', 'the value of MS1 tolerance';...
    'ms1_tolerance_type', 'the type of MS1 tolerance';...
    'ms2_tolerance', 'the value of MS2 tolerance';...
    'alpha', 'the relative intensity threshold using in denoising';...
    'fasta_file_path', 'the path of the fasta file';...
    'regular_express', 'the regular expression for parsing the protein names';...
    'model', 'the model used in algorithm (FEV/FEC/FEE)';...
    'method', 'the method used in algorithm (OLS/lasso)';...
    'result_filter_threshold', 'the result filter threshold of relative abundance';...
    'enzyme_name', 'the name of enzyme';...
    'enzyme_limit_C_term_possible_mod', 'the feasible C-terminal modification defined by the enzyme';...
    'output_dir_path', 'the path of output directory';...
    'checked_peptides_res_path', 'the path of checked peptide level result';...
    'msms_res_path', 'the path of msms level result'
};

% Check if all the required keys are present in the task_param_map
check_keys_in_map(task_param_map, key_needed_name_explain);
end



function check_required_keys_peptide_only(task_param_map)
% Check if all the required keys are present in the task_param_map for peptide only processing.
% Input:
%   task_param_map (containers.Map)
%       task parameter map

key_needed_name_explain = { ...
    % Parameter on peptide only processing
    'mod_file_path', 'the path of the modification database file';...
    'fixed_mod', 'the fixed modifications';...
    'variable_mod', 'the variable modifications';...
    'ms1_tolerance_value', 'the value of MS1 tolerance';...
    'ms1_tolerance_type', 'the type of MS1 tolerance';...
    'ms2_tolerance', 'the value of MS2 tolerance';...
    'alpha', 'the relative intensity threshold using in denoising';...
    'fasta_file_path', 'the path of the fasta file';...
    'regular_express', 'the regular expression for parsing the protein names';...
    'pep_spec_file_path', 'the path of the peptide-spectra list file';...
    'model', 'the model used in algorithm (FEV/FEC/FEE)';...
    'method', 'the method used in algorithm (OLS/lasso)';...
    'result_filter_threshold', 'the result filter threshold of relative abundance';...
    'enzyme_name', 'the name of enzyme';...
    'enzyme_limit_C_term_possible_mod', 'the feasible C-terminal modification defined by the enzyme';...
    'output_dir_path', 'the path of output directory';...
    'msms_res_path', 'the path of msms level result'
};

% Check if all the required keys are present in the task_param_map
check_keys_in_map(task_param_map, key_needed_name_explain);
end



function check_required_keys_site_level(task_param_map)
% Check if all the required keys are present in the task_param_map for site level processing.
% Input:
%   task_param_map (containers.Map)
%       task parameter map

key_needed_name_explain = { ...
    % Parameter on site level processing
    'protein_name_abbr_num','the number of protein name abbreviation';...
    'mod_name_abbr_num','the number of modification name abbreviation';...
    'pep_level_file_path','the path of peptide level result file';...
    'output_intere_path','the output path for the site-level result file of the target modification site';...
    'output_unintere_path','the output path for the peptide-level result file of the non-target modification site';...
    'ignore_strings_site_level','the strings that should be ignored in the peptide sequence to meet the need of isotopic labeled quantification methods'
};

% Check if all the required keys are present in the task_param_map
check_keys_in_map(task_param_map, key_needed_name_explain);

% 'protein_name_abbr_X' is the containers.Map of protein 'name' to 'abbreviation',
% where X is an integer between 1 and 'protein_name_abbr_num'.
% The following checks if all 'protein_name_abbr_X' are present.
for idx_pnan = 1:str2double(task_param_map('protein_name_abbr_num'))
    if ~task_param_map.isKey(['protein_name_abbr_',num2str(idx_pnan)])
        error(['The param ''protein_name_abbr_',num2str(idx_pnan),''' is not found.'])
    end
end

% 'mod_name_abbr_X' is the containers.Map of modification 'name' to 'abbreviation',
% where X is an integer between 1 and 'mod_name_abbr_num'.
% The following checks if all 'mod_name_abbr_X' are present.
for idx_mnan = 1:str2double(task_param_map('mod_name_abbr_num'))
    if ~task_param_map.isKey(['mod_name_abbr_',num2str(idx_mnan)])
        error(['The param ''mod_name_abbr_',num2str(idx_mnan),''' is not found.'])
    end
end
end



function check_required_keys_merge_to_pair_level(task_param_map)
% Check if all the required keys are present in the task_param_map for merge to pair level processing.
% Input:
%   task_param_map (containers.Map)
%       task parameter map

key_needed_name_explain = { ...
    % Parameter on merge to pair level processing
    'left_right_out_num','the number of pairs to be compared';...
    'left_name','the name of experiment on the left side';...
    'right_name','the name of experiment on the right side';...
    'ignore_strings_pair_level','the strings that should be ignored in the peptide sequence to meet the need of isotopic labeled quantification methods'
};

% Check if all the required keys are present in the task_param_map
check_keys_in_map(task_param_map, key_needed_name_explain);

% 'left_right_out_X' is a n*3 cell array, where n is the number of pairs to be compared,
% and each row contains the left experiment name, the right experiment name, and the output file name.
% The following checks if all 'left_right_out_X' are present.
for idx_lron = 1:str2double(task_param_map('left_right_out_num'))
    if ~task_param_map.isKey(['left_right_out_',num2str(idx_lron)])
        error(['The param ''left_right_out_',num2str(idx_lron),''' is not found.'])
    end
end
end



function check_required_keys_merge_pairs_level(task_param_map)
% Check if all the required keys are present in the task_param_map for merge pairs level processing.
% Input:
%   task_param_map (containers.Map)
%       task parameter map

key_needed_name_explain = { ...
    % Parameter on merge pairs level processing
    'pair_num','the number of pairs to be merged';...
    'final_output_path','the output path for the merged result file'
};

% Check if all the required keys are present in the task_param_map
check_keys_in_map(task_param_map, key_needed_name_explain);

% 'pair_X' is the path of the result file of the X-th pair to be merged.
% The following checks if all 'pair_X' are present.
for idx_pair = 1:str2double(task_param_map('pair_num'))
    if ~task_param_map.isKey(['pair_',num2str(idx_pair)])
        error(['The param ''pair_',num2str(idx_pair),''' is not found.'])
    end
end

% 'left_right_name_X' is a n*2 cell array indicating the name of the X-th pair to be merged.
% The following checks if all 'left_right_name_X' are present.
for idx_pair = 1:str2double(task_param_map('pair_num'))
    if ~task_param_map.isKey(['left_right_name_',num2str(idx_pair)])
        error(['The param ''left_right_name_',num2str(idx_pair),''' is not found.'])
    end
end
end



function check_keys_in_map(task_param_map, key_needed_name_explain)
% Check if all the required keys are present in the task_param_map.
% Input:
%   task_param_map
%       A containers.Map object that contains the task parameters.
%   key_needed_name_explain
%       A cell array where each row contains the key name and the explanation of the key.

for iKey = 1:size(key_needed_name_explain,1)
    if ~task_param_map.isKey(key_needed_name_explain{iKey, 1})
        error(['The param ''',key_needed_name_explain{iKey, 1},''' which specifies ', ...
            key_needed_name_explain{iKey, 2},' is not found.'])
    end
end
end