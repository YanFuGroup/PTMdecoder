function res_joint = get_merge_res_of_pairs(obj)
% Get the merge result of pairs
% Output:
%   res_joint
%       Merged result of pairs in site-level with joint of peptides in all pairs.
%           For example, {H3K4ac, 
%                           ['PEPT{Acetyl}IDE','+2',{},{1e5,2e5}; 'PEPT{Acetyl}IDE','+3',{3e6,4e6},{}] }

% Read the merge result of pairs
res_content = cell(length(obj.m_result_paths), 1);
for i = 1:length(obj.m_result_paths)
    res_content{i} = read_each_pair(obj.m_result_paths{i}, obj.m_column_idxs);
end

% Merge the results of pairs
res_joint = get_joint_res_content(res_content);
end



%% Other functions
function res_content = read_each_pair(result_path, column_idxs)
% Read the merge result of a pair
% Input:
%   result_path
%       Path of the merge result of a pair
%   column_idxs
%       Indexes of columns to read
% Output:
%   res_content
%       Content of the merge result of a pair

% Open the file
fin = fopen(result_path, 'r');
if fin == -1
    error('Cannot open the file of pair: "%s"!', result_path);
end

% Read the file and check how many sites are there
site_nums = 0;
strline = fgetl(fin);   %#ok<NASGU> % Skip the header
while ~feof(fin)
    strline = fgetl(fin);
    if ~isequal(strline(1), char(9))    % Not '\t', is a new site
        site_nums = site_nums + 1;
    end
end

% Rewind the file
frewind(fin);

% Initialize the result [site, ModifiedPeptideSet]
res_content = cell(site_nums, 2);

% Read the file and store the content
site_idx = 0;
strline = fgetl(fin);   %#ok<NASGU> % Skip the header
while ~feof(fin)
    strline = fgetl(fin);
    if isempty(strline)
        continue;
    end
    segment = strsplit(strline, '\t');
    if ~isequal(strline(1), char(9))    % Not '\t', is a new site
        % A new site
        site_idx = site_idx + 1;
        res_content{site_idx, 1} = segment{column_idxs.icol_site};
        res_content{site_idx, 2} = {segment{column_idxs.icol_pep}, segment{column_idxs.icol_charge}, ...
            segment{column_idxs.icol_quant_1}, segment{column_idxs.icol_quant_2}};
    else
        % Append the content of the last site
        res_content{site_idx, 2} = [res_content{site_idx, 2}; ...
            {segment{column_idxs.icol_pep}, segment{column_idxs.icol_charge}, ...
            segment{column_idxs.icol_quant_1}, segment{column_idxs.icol_quant_2}}];
    end
end

% Close the file
fclose(fin);
end



function res_joint = get_joint_res_content(res_content)
% Merge the results of pairs
% Input:
%   res_content
%       Content of the merge result of pairs
% Output:
%   res_joint
%       Merge result of pairs, joint of peptides in all pairs

% Initialize the result, [site, ModifiedPeptideSet]
site_joint = cellfun(@(x) [x(:,1)], res_content, 'UniformOutput', false);
site_joint = unique(vertcat(site_joint{:})); % Unique the joint sites
res_joint = [site_joint, cell(size(site_joint, 1), 1)];

% Transverse the sites and merge the peptides
for idx_site = 1:size(site_joint, 1)
    % Extract the peptides of current site
    extracted_peptides = get_peptides_of_site(res_content, site_joint{idx_site, 1});

    % Merge the peptides
    merged_peptides = merge_peptides(extracted_peptides);

    % Store the merged peptides
    res_joint{idx_site, 2} = merged_peptides;
end
end



function extracted_peptides = get_peptides_of_site(res_content, site_name)
% Extract the peptides of a site from res_content
% Input:
%   res_content
%       Content of the merge result of pairs
%   site_name
%       Site to extract
% Output:
%   extracted_peptides
%       Peptides of the site, cell array. Each element is a pair of experiment. 
%       Each pair is a cell array, each element is a peptide. Empty if no peptide.

% Initialize the result
extracted_peptides = cell(length(res_content), 1);

% Transverse the pairs and extract the peptides
for idx_pair = 1:length(res_content)
    % Extract the peptides of current site
    idx_site = find(ismember(res_content{idx_pair}(:,1), site_name));
    if ~isempty(idx_site)
        extracted_peptides{idx_pair} = res_content{idx_pair}{idx_site, 2};
    end
end
end



function merged_peptides = merge_peptides(extracted_peptides)
% Merge the peptides of a site
% Input:
%   extracted_peptides
%       Peptides of the site, cell array. Each element is a pair of experiment.
%       Each pair is a cell array, each element is a peptide. Empty if no peptide.
% Output:
%   merged_peptides
%       Merged peptides, cell array. Each element is a pair of quantification as numeric. Empty if no peptide exist in the pair.
%           For example, {'PEPTIDE','+2',[],[1e5,2e5]; 'PEPTIDE','+3',[3e6,4e6],[]}

% Concatenate the first column (peptide sequence) and the second column (charge) 
%   in each extracted_peptides to form a 'peptide sequence + charge' cell array, then remove duplicates.
pep_charge_cell = cell(length(extracted_peptides), 1);
for idx_pair = 1:length(extracted_peptides)
    if isempty(extracted_peptides{idx_pair}) || size(extracted_peptides{idx_pair}, 2) < 2
        pep_charge_cell{idx_pair} = [];
    else
        pep_charge_cell{idx_pair} = extracted_peptides{idx_pair}(:, 1:2);
    end
end
reshaped_pep_charge = vertcat(pep_charge_cell{:});
pep_charge_joint = strcat(reshaped_pep_charge(:,1),reshaped_pep_charge(:,2));
[pep_charge_joint, uni_idxs] = unique(pep_charge_joint); % Unique the joint peptide-charge chars
uni_pep_charge = reshaped_pep_charge(uni_idxs, :);

% Initialize the result
merged_peptides = [uni_pep_charge, cell(size(pep_charge_joint, 1), size(extracted_peptides, 1))];

% Transverse the joint peptide-charge chars and merge the quantifications of the same peptide
for idx_pep = 1:length(pep_charge_joint)
    for idx_pair = 1:length(extracted_peptides)
        if isempty(extracted_peptides{idx_pair})
            continue;
        end
        idx_pep_in_pair = find(ismember(strcat(extracted_peptides{idx_pair}(:,1), ...
            extracted_peptides{idx_pair}(:,2)), pep_charge_joint{idx_pep, 1}));
        if ~isempty(idx_pep_in_pair)
            % Merge the quantifications of the same peptide
            merged_peptides{idx_pep, 2 + idx_pair} = extracted_peptides{idx_pair}(idx_pep_in_pair, 3:4);
        end
    end
end
end