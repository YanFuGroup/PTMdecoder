function result = FilterWithChemPrior(result)
% Filter the results before FDR control
% The chemical prior: 
%   - There will only be at most two modification types on a peptide. 
%   - Only some modifications are allowed at C-terminus of a peptide.
% input & output:
%   result
%       the result read from the search engine result file


%% Settings of prior knowledge
max_mod_num = 2; % the maximum number of modification types on a peptide
possible_C_term_mod = {'Methyl', 'Label:13C(6)', 'Label:13C(6)15N(4)', 'Label:13(6)+Methyl', 'Label:13(6)15(4)+Methyl'}; % the allowed modifications at C-terminus of a peptide


%% Filter the results

% Initialize a logical array to store filtering results
filter_mask = true(size(result));

for idx = 1:length(result)
    mod_name = result(idx).modification;
    mod_pos = result(idx).modificationlocation;
    pep = result(idx).peptide;

    is_del = false;
    if ~isequal(mod_name, '-')
        %% Check for at most two modification types
        str_cell_temp = strsplit(mod_name, ',');
        str_cell_temp = cellfun(@(y) regexp(y, '^[^\s]+', 'match'), str_cell_temp, 'UniformOutput', false);
        unique_mod = unique([str_cell_temp{:}]); % Extract unique modifications
        if length(unique_mod) > max_mod_num
            is_del = true; % Skip this result if it exceeds the max modification number
        end

        %% Only some modifications are allowed at C-terminus of a peptide.
        mod_positions = strsplit(mod_pos, ',');
        max_val = max(str2double(mod_positions));
        % Check if the maximum modification position is equal to the length of the peptide
        if max_val == length(pep)
            % Extract the name of the C-terminal (amino acid) modification
            mod_positions_num = str2double(mod_positions);
            max_indices = find(mod_positions_num == max_val);
            mod_names = strsplit(mod_name, ',');
            for max_idx = max_indices
                mod_names{max_idx} = strsplit(mod_names{max_idx}, ' ');
                mod_names{max_idx} = mod_names{max_idx}{1};

                % Check if the modification name is not in the allowed list
                if ~any(ismember(mod_names{max_idx}, possible_C_term_mod))
                    is_del = true;
                end
            end
        end
    end

    % Update the filter mask
    if is_del
        filter_mask(idx) = false;
    end
end

% Apply the filter mask to the result
result = result(filter_mask);

end

