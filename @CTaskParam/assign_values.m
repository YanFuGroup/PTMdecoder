function obj = assign_values(obj, task_param_map)
% This function assigns values to the properties of the CTaskParam object based on the provided task_param_map.
% Input:
%   task_param_map
%       A map containing the parameter names and their corresponding values.

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
for iKey = 1:length(module_names)
    if task_param_map.isKey(module_names{iKey}) && strcmp(task_param_map(module_names{iKey}), '1')
        obj.(['m_', module_names{iKey}]) = true;
    else
        obj.(['m_', module_names{iKey}]) = false;
    end
end

% MSMS and peptide level
if obj.m_msms_peptide_level_on
    obj = set_msms_peptide_level_params_from_tpmap(obj, task_param_map);

% Peptide requantification
elseif obj.m_peptide_requant_on
    obj = set_peptide_requant_params_from_tpmap(obj, task_param_map);

% Peptide only
elseif obj.m_peptide_only_on
    obj = set_peptide_only_params_from_tpmap(obj, task_param_map);
end

% The following processes can be run independently

% Site level
if obj.m_site_level_on
    obj = set_site_level_params_from_tpmap(obj, task_param_map);
end

% Merge to pair level
if obj.m_merge_to_pair_level_on
    obj = set_merge_to_pair_level_params_from_tpmap(obj, task_param_map);
end

% Merge pairs level
if obj.m_merge_pairs_level_on
    obj = set_merge_pairs_level_params_from_tpmap(obj, task_param_map);
end
end



%% Other functions

function obj = set_msms_peptide_level_params_from_tpmap(obj, task_param_map)
% This function assigns values to the properties of the CTaskParam object based on the provided task_param_map.
% Input:
%   task_param_map
%       A map containing the parameter names and their corresponding values.

% About modifications, including fixed_mod, variable_mod, and mod_file
obj.m_mod_file_path = task_param_map('mod_file_path');
obj.m_fixed_mod = task_param_map('fixed_mod');
obj.m_variable_mod = task_param_map('variable_mod');

% About the spectra file, including the path, tolerance, and filter threshold
obj.m_spec_dir_path = task_param_map('spec_dir_path');
obj.m_ms1_tolerance.value = str2double(task_param_map('ms1_tolerance_value'));
% Check if ms1_tolerance_type is 'PPM' and set isppm accordingly
if isequal(upper(task_param_map('ms1_tolerance_type')),'PPM')
    obj.m_ms1_tolerance.isppm = 1;
else
    obj.m_ms1_tolerance.isppm = 0;
end
obj.m_ms2_tolerance = str2double(task_param_map('ms2_tolerance'));
obj.m_alpha = str2double(task_param_map('alpha'));

% About the database file, including the path and regular expression
obj.m_fasta_file_path = task_param_map('fasta_file_path');
obj.m_regular_express = task_param_map('regular_express');

% About the pep_spec file, including the path
obj.m_pep_spec_file_path = task_param_map('pep_spec_file_path');

% About the parameters for the core algorithm
obj.m_model = str2double(task_param_map('model'));
obj.m_method = str2double(task_param_map('method'));
obj.m_lambda = str2double(task_param_map('lambda'));
obj.m_result_filter_threshold = str2double(task_param_map('result_filter_threshold'));

% About the enzyme, including the name and limits, for specific experiments
obj.m_enzyme_name = task_param_map('enzyme_name');
obj.m_enzyme_limits = str2num(task_param_map('enzyme_limit_C_term_possible_mod')); %#ok<ST2NM> 

% About the output file, including the path
obj.m_output_dir_path = task_param_map('output_dir_path');

% About min_MSMS_num parameter
if task_param_map.isKey('min_MSMS_num')
    obj.m_min_MSMS_num = str2double(task_param_map('min_MSMS_num'));
else
    obj.m_min_MSMS_num = 1; % Default value
end

% About the checked peptides result path and msms result path (not used in this process)
obj.m_checked_peptides_res_path = [];
obj.m_msms_res_path = [];
end



function obj = set_peptide_requant_params_from_tpmap(obj, task_param_map)
% This function assigns values to the properties of the CTaskParam object based on the provided task_param_map.
% Input:
%   task_param_map
%       A map containing the parameter names and their corresponding values.

% About modifications, including fixed_mod, variable_mod, and mod_file
obj.m_mod_file_path = task_param_map('mod_file_path');
obj.m_fixed_mod = task_param_map('fixed_mod');
obj.m_variable_mod = task_param_map('variable_mod');

% About the spectra file, including the path, tolerance, and filter threshold
obj.m_spec_dir_path = task_param_map('spec_dir_path');
obj.m_ms1_tolerance.value = str2double(task_param_map('ms1_tolerance_value'));
% Check if ms1_tolerance_type is 'PPM' and set isppm accordingly
if isequal(upper(task_param_map('ms1_tolerance_type')),'PPM')
    obj.m_ms1_tolerance.isppm = 1;
else
    obj.m_ms1_tolerance.isppm = 0;
end
obj.m_ms2_tolerance = str2double(task_param_map('ms2_tolerance'));
obj.m_alpha = str2double(task_param_map('alpha'));

% About the database file, including the path and regular expression
obj.m_fasta_file_path = task_param_map('fasta_file_path');
obj.m_regular_express = task_param_map('regular_express');

% About the pep_spec file, including the path
% obj.m_pep_spec_file_path = task_param_map('pep_spec_file_path');
obj.m_pep_spec_file_path = [];

% About the parameters for the core algorithm
obj.m_model = str2double(task_param_map('model'));
obj.m_method = str2double(task_param_map('method'));
obj.m_lambda = str2double(task_param_map('lambda'));
obj.m_result_filter_threshold = str2double(task_param_map('result_filter_threshold'));

% About the enzyme, including the name and limits, for specific experiments
obj.m_enzyme_name = task_param_map('enzyme_name');
obj.m_enzyme_limits = str2num(task_param_map('enzyme_limit_C_term_possible_mod')); %#ok<ST2NM> 

% About the output file, including the path
obj.m_output_dir_path = task_param_map('output_dir_path');

% About min_MSMS_num parameter
if task_param_map.isKey('min_MSMS_num')
    obj.m_min_MSMS_num = str2double(task_param_map('min_MSMS_num'));
else
    obj.m_min_MSMS_num = 1; % Default value
end

% About the checked peptides result path and msms result path
obj.m_checked_peptides_res_path = task_param_map('checked_peptides_res_path');
obj.m_msms_res_path = task_param_map('msms_res_path');
end



function obj = set_peptide_only_params_from_tpmap(obj, task_param_map)
% This function assigns values to the properties of the CTaskParam object based on the provided task_param_map.
% Input:
%   task_param_map
%       A map containing the parameter names and their corresponding values.

% About modifications, including fixed_mod, variable_mod, and mod_file
obj.m_mod_file_path = task_param_map('mod_file_path');
obj.m_fixed_mod = task_param_map('fixed_mod');
obj.m_variable_mod = task_param_map('variable_mod');

% About the spectra file, including the path, tolerance, and filter threshold
obj.m_spec_dir_path = task_param_map('spec_dir_path');
obj.m_ms1_tolerance.value = str2double(task_param_map('ms1_tolerance_value'));
% Check if ms1_tolerance_type is 'PPM' and set isppm accordingly
if isequal(upper(task_param_map('ms1_tolerance_type')),'PPM')
    obj.m_ms1_tolerance.isppm = 1;
else
    obj.m_ms1_tolerance.isppm = 0;
end
obj.m_ms2_tolerance = str2double(task_param_map('ms2_tolerance'));
obj.m_alpha = str2double(task_param_map('alpha'));

% About the database file, including the path and regular expression
obj.m_fasta_file_path = task_param_map('fasta_file_path');
obj.m_regular_express = task_param_map('regular_express');

% About the pep_spec file, including the path
obj.m_pep_spec_file_path = [];

% About the parameters for the core algorithm
obj.m_model = str2double(task_param_map('model'));
obj.m_method = str2double(task_param_map('method'));
obj.m_lambda = str2double(task_param_map('lambda'));
obj.m_result_filter_threshold = str2double(task_param_map('result_filter_threshold'));

% About the enzyme, including the name and limits, for specific experiments
obj.m_enzyme_name = task_param_map('enzyme_name');
obj.m_enzyme_limits = str2num(task_param_map('enzyme_limit_C_term_possible_mod')); %#ok<ST2NM> 

% About the output file, including the path
obj.m_output_dir_path = task_param_map('output_dir_path');

% About the checked peptides result path and msms result path
obj.m_checked_peptides_res_path = [];
obj.m_msms_res_path = task_param_map('msms_res_path');
end



function obj = set_site_level_params_from_tpmap(obj, task_param_map)
% This function assigns values to the properties of the CTaskParam object based on the provided task_param_map.
% Input:
%   task_param_map
%       A map containing the parameter names and their corresponding values.

% About the target protein name and their abbreviation
obj.m_protein_name_abbr_num = str2double(task_param_map('protein_name_abbr_num'));
obj.m_protein_name_abbr = containers.Map;
for idx_pnan = 1:obj.m_protein_name_abbr_num
    protein_name_abbr_str = task_param_map(['protein_name_abbr_',num2str(idx_pnan)]);
    split_str = strsplit(protein_name_abbr_str, '>');
    protein_name = strtrim(split_str{1});
    protein_abbr = strtrim(split_str{2});
    obj.m_protein_name_abbr(protein_name) = protein_abbr;
end

% About the target modification name and their abbreviation
obj.m_mod_name_abbr_num = str2double(task_param_map('mod_name_abbr_num'));
obj.m_mod_name_abbr = containers.Map;
for idx_mnan = 1:obj.m_mod_name_abbr_num
    mod_name_abbr_str = task_param_map(['mod_name_abbr_',num2str(idx_mnan)]);
    split_str = strsplit(mod_name_abbr_str, '>');
    mod_name = strtrim(split_str{1});
    mod_abbr = strtrim(split_str{2});
    obj.m_mod_name_abbr(mod_name) = mod_abbr;
end

% About the pathes
three_paths_name = {'pep_level_file_path', 'output_intere_path', 'output_unintere_path'};
three_paths_default = {'report_peptide_all.txt', 'report_site.txt', 'report_peptide_uninterested.txt'};
for i_tpn = 1:length(three_paths_name)
    if ~isempty(task_param_map(three_paths_name{i_tpn}))
        obj.(['m_',three_paths_name{i_tpn}]) = task_param_map(three_paths_name{i_tpn});
    else
        if ~isempty(task_param_map('output_dir_path'))
            obj.(['m_',three_paths_name{i_tpn}]) = fullfile(task_param_map('output_dir_path'), three_paths_default{i_tpn});
        end
    end
end

% Ignore the following strings in the modified peptide sequence strings,
% to meet the need of analysis of labeled quantification methods such as SILAC
obj.m_ignore_strings_site_level = get_ignore_strings(task_param_map('ignore_strings_site_level'));
end



function obj = set_merge_to_pair_level_params_from_tpmap(obj, task_param_map)
% This function assigns values to the properties of the CTaskParam object based on the provided task_param_map.
% Input:
%   task_param_map
%       A map containing the parameter names and their corresponding values.

% About the pathes
obj.m_left_right_out_num = str2double(task_param_map('left_right_out_num'));
obj.m_left_right_out = cell(obj.m_left_right_out_num, 3);
for idx_lron = 1:obj.m_left_right_out_num
    % The format of the string is 'left|right>out'
    left_right_out_str = task_param_map(['left_right_out_',num2str(idx_lron)]);
    split_str = strsplit(left_right_out_str, {'|', '>'});
    obj.m_left_right_out(idx_lron, :) = strtrim(split_str);
end

% About the names
obj.m_left_name = task_param_map('left_name');
obj.m_right_name = task_param_map('right_name');

% Ignore the following strings in the modified peptide sequence strings,
% to meet the need of analysis of labeled quantification methods such as SILAC
obj.m_ignore_strings_pair_level = get_ignore_strings(task_param_map('ignore_strings_pair_level'));
end



function obj = set_merge_pairs_level_params_from_tpmap(obj, task_param_map)
% This function assigns values to the properties of the CTaskParam object based on the provided task_param_map.
% Input:
%   task_param_map
%       A map containing the parameter names and their corresponding values.

% About the pathes and names of each pair
obj.m_pair_num = str2double(task_param_map('pair_num'));
obj.m_pair = cell(obj.m_pair_num, 1);
obj.m_left_right_name = cell(obj.m_pair_num, 2);
for idx_pn = 1:obj.m_pair_num
    % The format of the 'left_right_name' is 'left|right'
    obj.m_pair{idx_pn} = task_param_map(['pair_',num2str(idx_pn)]);
    left_right_name_str = task_param_map(['left_right_name_',num2str(idx_pn)]);
    split_str = strsplit(left_right_name_str, '|');
    obj.m_left_right_name{idx_pn, 1} = strtrim(split_str{1});
    obj.m_left_right_name{idx_pn, 2} = strtrim(split_str{2});
end

% Path of the output file
obj.m_final_output_path = task_param_map('final_output_path');
end



function ignore_strings = get_ignore_strings(ignore_strings_str)
% This function converts the ignore_strings_str to a cell array of char strings.
% Input:
%   ignore_strings_str
%       A string containing the ignore strings, separated by ';'.
%       Each ignore string is enclosed in double quotes.
%       Example: '"string1";"string2";"string3"'
% Output:
%   ignore_strings
%       A cell array of char strings containing the ignore strings.

ignore_strings = [];
ignore_strings_seg = strsplit(ignore_strings_str, ';');

% Remove the double quotes from the ignore strings
ignore_strings_seg = cellfun(@(s) regexp(s, '"([^"]*)"', 'tokens'), ignore_strings_seg, 'UniformOutput', false);
for idx_is = 1:length(ignore_strings_seg)
    if isempty(ignore_strings_seg{1})
        continue;
    else
        ignore_strings = [ignore_strings, ignore_strings_seg{idx_is}{1}]; %#ok<AGROW>
    end
end
end