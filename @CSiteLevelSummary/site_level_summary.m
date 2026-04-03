function obj = site_level_summary(obj)
% Summary the quantification of each modified PSM on site level.
% Input:
%   obj (CSiteLevelSummary)
%       site-level summarizer instance
% Using members of obj:
%   m_input_path
%       the path of the result file on peptide level
%   m_protein_name_abbr
%       the map of target protein name and abbreviation
%   m_mod_name_abbr
%       the map of target modification name and abbreviation
%   m_ignore_strings
%       the strings to be ignored in peptide sequence string, 
%       skip the fixed modification
%   m_column_idxs
%       the index of the columns in the result file
% Update members of obj:
%   m_result_output_index
%       the map of result site, such as "H1K1ac", and the index of it in
%       the following two cell
%   m_result_output_string
%       the cell of output string corresbonding to the target site, the 
%       index is record in result_output_index
%   m_result_output_sum
%       the cell of sum of intensity corresbonding to the target site, the
%       index is record in result_output_index
%   m_result_uninterested_string
%       the cell of uninterested corresbonding to the target site, the
%       index is record in result_output_index

% site-level summary map, for searching the string
% [site string] -> [index number of the output string cell]
result_output_index = containers.Map;
buff_length = 1000;
result_output_string = cell(buff_length, 1);
result_output_sum = cell(buff_length, 1);
result_output_length = 0;
result_uninterested_string = cell(buff_length, 1);
result_uninterested_length = 0;

total_protein_lines = 0;
matched_protein_lines = 0;
unmatched_protein_lines = 0;
uninterested_peptide_lines = 0;

% open files and check
fin = fopen(obj.m_input_path);
if fin < 0
    error(['Can not open the peptide level result file: "', obj.m_input_path, '"']);
end

CLogger.info('[CSiteLevelSummary:site_level_summary] Start summarizing. input=%s', obj.m_input_path);

% read files and gather the specified sequence
fgetl(fin);
selected_abbr_protein = [];
fgetl(fin); % skip the first two lines, they are header lines
while ~feof(fin)
    strline = fgetl(fin);
    if isempty(strline)
        continue;
    end

    % read one line and process
    if ~isequal(strline(1),'*') && ~isequal(strline(1),'@')
        % the protein lines
        total_protein_lines = total_protein_lines + 1;
        selected_abbr_protein = [];
        selected_start_pos_protein = -1;
        segments = split(strline, ';');
        for i = 1:length(segments)-1
            % Split the part into key and number
            keyValue = split(segments{i}, ',');
            % Extract the key and convert it to a string
            if isKey(obj.m_protein_name_abbr, char(keyValue(1)))
                % Keep legacy strategy: use the first resolvable hit.
                selected_abbr_protein = obj.m_protein_name_abbr(char(keyValue(1)));
                selected_start_pos_protein = str2double(char(keyValue(2)));
                break;
            end
        end
        % no target protein, record as the uninterested protein
        if isempty(selected_abbr_protein)
            unmatched_protein_lines = unmatched_protein_lines + 1;
            % record the uninterested string
            result_uninterested_length = result_uninterested_length + 1;
            if result_uninterested_length > length(result_uninterested_string)
                result_uninterested_string{result_uninterested_length+buff_length} = [];
            end
            result_uninterested_string{result_uninterested_length} = strline;
        else
            matched_protein_lines = matched_protein_lines + 1;
        end
    elseif isequal(strline(1), '*')
        % check if the protein is found
        if ~isempty(selected_abbr_protein)
            % the IMP and quantification lines
            segments = split(strline);
            
            % process the string of IMP
            modified_peptides = segments{obj.m_column_idxs.icol_seq};

            % ignore the obj.m_ignore_strings
            for idx_ig = 1:length(obj.m_ignore_strings)
                modified_peptides = strrep(modified_peptides, obj.m_ignore_strings{idx_ig}, '');
            end

            % Find all substrings between curly braces
            mod_str_matches = regexp(modified_peptides, '{(.*?)}', 'tokens');

            % Extract the substrings and positions
            mods_str = cell(1, numel(mod_str_matches));
            positions_seq = zeros(1, numel(mod_str_matches));
            positions_str = zeros(1, numel(mod_str_matches));
            mod_specificity = cell(1, numel(mod_str_matches));
            mod_prot_pos = repmat(-1, 1, numel(mod_str_matches));
            start_pos = 0;
            for i = 1:numel(mod_str_matches)
                % matched modification string
                mod_name = mod_str_matches{i}{1};
                found_index = strfind(modified_peptides(start_pos+1:end), ['{' mod_name '}']);
                if i == 1
                    positions_seq(i) = found_index(1) - 2;
                else
                    positions_seq(i) = positions_seq(i-1) + found_index(1) - 1;
                end
                % accumulate the position and ready to find next
                positions_str(i) = start_pos + found_index(1);
                mods_str{i} = mod_name;
                start_pos = positions_str(i) + numel(mod_name) + 1;

                % whether in the map, whether need to be collect
                if isKey(obj.m_mod_name_abbr, mod_name)
                    abbr_mod = obj.m_mod_name_abbr(mod_name);
                else
%                     abbr_mod = mod_name;
                    % if the modification is not found in the
                    % 'obj.m_mod_name_abbr' map, then do not record it
                    continue;
                end

                % record the specificity of the modification
                mod_specificity{i} = modified_peptides(positions_str(i)-1);
                if mod_specificity{i} == '_'
                    % change '_' to corresbonding N/C-term
                    if positions_seq(i) == 0
                        mod_specificity{i} = 'N-term';
                    else
                        mod_specificity{i} = 'C-term';
                    end
                    site_name = [selected_abbr_protein, mod_specificity{i}, '_', abbr_mod];
                else
                    % calculate the position of site on the protein
                    if selected_start_pos_protein < 0
                        error(['The start position on protein of peptide "', ...
                            modified_peptides, '" is out of range.']);
                    end
                    mod_prot_pos(i) = selected_start_pos_protein + positions_seq(i) - 1;
                    if obj.m_site_position_count_initial_m
                        site_pos = mod_prot_pos(i);
                    else
                        site_pos = mod_prot_pos(i) - 1;
                    end
                    site_name = [selected_abbr_protein, mod_specificity{i} , num2str(site_pos), abbr_mod];
                end

                % append the line to the corresbonding result set
                if isKey(result_output_index, site_name)
                    % append to the existing result cell
                    result_output_string{result_output_index(site_name)} = ...
                        [result_output_string{result_output_index(site_name)}, strline, '\n'];
                    result_output_sum{result_output_index(site_name)} = ...
                        result_output_sum{result_output_index(site_name)} + str2double(segments{obj.m_column_idxs.icol_auc});
                else
                    % add a new site result map key
                    result_output_length = result_output_length + 1;
                    result_output_index(site_name) = result_output_length;
                    if result_output_length > length(result_output_string)
                        % if the capacity is not enough, alloc some space
                        result_output_string{result_output_length+buff_length} = [];
                        result_output_sum{result_output_length+buff_length} = [];
                    end
                    result_output_string{result_output_length} = [strline, '\n'];
                    result_output_sum{result_output_length} = str2double(segments{obj.m_column_idxs.icol_auc});
                end
            end
        else
            % if the protein is not interested, then record it in cell
            uninterested_peptide_lines = uninterested_peptide_lines + 1;
            result_uninterested_length = result_uninterested_length + 1;
            if result_uninterested_length > length(result_uninterested_string)
                result_uninterested_string{result_uninterested_length+buff_length} = [];
            end
            result_uninterested_string{result_uninterested_length} = strline;
        end
    end
end
result_output_string(result_output_length+1:end) = [];
result_output_sum(result_output_length+1:end) = [];
result_uninterested_string(result_uninterested_length+1:end) = [];

fclose(fin);

% update the members of obj
obj.m_result_output_index = result_output_index;
obj.m_result_output_string = result_output_string;
obj.m_result_output_sum = result_output_sum;
obj.m_result_uninterested_string = result_uninterested_string;

if unmatched_protein_lines > 0
    CLogger.warn(['[CSiteLevelSummary:site_level_summary] Unmatched protein lines found: %d ', ...
        '(mapped=%d, total=%d).'], unmatched_protein_lines, matched_protein_lines, total_protein_lines);
end

if uninterested_peptide_lines > 0
    CLogger.warn('[CSiteLevelSummary:site_level_summary] Uninterested peptide lines recorded: %d.', uninterested_peptide_lines);
end

CLogger.info(['[CSiteLevelSummary:site_level_summary] Done. total_protein_lines=%d, ', ...
    'matched_protein_lines=%d, interested_sites=%d, uninterested_records=%d.'], ...
    total_protein_lines, matched_protein_lines, result_output_length, result_uninterested_length);

end