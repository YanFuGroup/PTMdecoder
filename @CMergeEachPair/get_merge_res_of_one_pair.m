function res_intersec = get_merge_res_of_one_pair(obj)
% Merge the site-level result for each pair
% Output:
%   res_intersec
%       the merged result of the two samples

% Read the two files and record them in the format of [site_name, quantification, peptides(cell)]
res_content1 = read_site_file(obj.m_result_path1, obj.m_column_idxs);
res_content2 = read_site_file(obj.m_result_path2, obj.m_column_idxs);

% Compare the two result content and return the output
res_intersec = get_intersection_res_content(res_content1, res_content2, obj.m_column_idxs, obj.m_ignore_strings);
end



% Other functions

function res_intersec = get_intersection_res_content(res_content1, res_content2, column_idxs, ignore_strings)
% Get the intersection of two site-level result
% Input:
%   res_content1, res_content2
%       the site-level result of two samples
%   column_idxs
%       the index of columns of the site-level result
%       [icol_site, icol_pep, icol_charge, icol_dataset, icol_quant, icol_max]
%   ignore_strings
%       the titles to be ignored in the peptide sequence
% Output:
%   res_intersec
%       the intersection of the two site-level result

% Sort the two result content by site name
res_content1 = sortrows(res_content1,1);
res_content2 = sortrows(res_content2,1);

% Initialize the output
buff_length = 100;
res_intersec = cell(buff_length,2); % [site_name, intersection_result(sequnce, charge, quant1, quant2)]
idx1 = 1;
idx2 = 1;

% Compare the two result content and record the intersection
idx_intersec = 0;
while idx1 <= size(res_content1,1) && idx2 <= size(res_content2,1)
    res = str_cmp_dict(res_content1{idx1,1}, res_content2{idx2,1});
    if res == 0
        % The site name is equal, record the intersection
        idx_intersec = idx_intersec + 1;

        % If the buffer is not enough, expand the buffer
        if idx_intersec > size(res_intersec,1)
            res_intersec{idx_intersec+buff_length,1} = '';
        end

        % Record the intersection site
        res_intersec{idx_intersec,1} = res_content1{idx1,column_idxs.icol_site};
        res_intersec{idx_intersec,2} = get_intersection_peptides(res_content1{idx1,3}, ...
            res_content2{idx2,3}, column_idxs, ignore_strings);
        %         res_intersec{idx_intersec,2} = ;

        % If the intersection is empty, remove the row
        if isempty(res_intersec{idx_intersec,2})
            res_intersec{idx_intersec,1} = '';
            idx_intersec = idx_intersec - 1;
        end

        idx1 = idx1 + 1;
        idx2 = idx2 + 1;
    elseif res < 0
        % The site name of the first is smaller, move the first index
        idx1 = idx1 + 1;
    else
        % The site name of the second is smaller, move the second index
        idx2 = idx2 + 1;
    end
end

% Remove the empty rows
res_intersec(idx_intersec+1:end,:) = [];
end



function res_intersec_pep = get_intersection_peptides(pep1, pep2, column_idxs, ignore_strings)
% Get the intersection of two peptide list
% Input:
%   pep1, pep2
%       the peptide list of two samples
%   column_idxs
%       the index of columns of the peptide list
%       [icol_site, icol_pep, icol_charge, icol_dataset, icol_quant, icol_max]
%   ignore_strings
%       the titles to be ignored in the peptide sequence
% Output:
%   res_intersec_pep
%       the intersection of the two peptide list

% Initialize the output
res_intersec_pep = [];

% Organize the two peptide list, [name, sequence, charge, quantification]
pep1 = organize_peptide_list(pep1, column_idxs, ignore_strings);
pep2 = organize_peptide_list(pep2, column_idxs, ignore_strings);

% Sort the two peptide list by peptide sequence
pep1 = sortrows(pep1,1);
pep2 = sortrows(pep2,1);

% Compare the two peptide list and record the intersection
idx1 = 1;
idx2 = 1;
while idx1 <= size(pep1,1) && idx2 <= size(pep2,1)
    res_cmp = str_cmp_dict(pep1{idx1,1}, pep2{idx2,1});
    if res_cmp == 0
        % The peptide sequence is equal, record the intersection
        res = {pep1{idx1,2},pep1{idx1,3},pep1{idx1,4},pep2{idx2,4}};
        res_intersec_pep = [res_intersec_pep; res]; %#ok<AGROW>
        idx1 = idx1 + 1;
        idx2 = idx2 + 1;
    elseif res_cmp < 0
        % The peptide sequence of the first is smaller, move the first index
        idx1 = idx1 + 1;
    else
        % The peptide sequence of the second is smaller, move the second index
        idx2 = idx2 + 1;
    end
end
end



function res = organize_peptide_list(pep, column_idxs, ignore_strings)
% Organize the peptide list
% Input:
%   pep
%       the peptide list
%   column_idxs
%       the index of columns of the peptide list
%       [icol_site, icol_pep, icol_charge, icol_dataset, icol_quant, icol_max]
%   ignore_strings
%       the titles to be ignored in the peptide sequence
% Output:
%   res
%       the organized peptide list, [name, sequence, charge, quantification]

% Initialize the output, [name, sequence, charge, quantification]
res = cell(size(pep,1),4);

% Record the peptide list
for i = 1:size(pep,1)
    % Split the peptide line
    segment = split(pep{i});

    % Delete the ignored strings
    temp_pep = segment{column_idxs.icol_pep};
    for idx_is = 1:length(ignore_strings)
        temp_pep = strrep(temp_pep, ignore_strings{idx_is}, '');
    end
    res{i,1} = [temp_pep,segment{column_idxs.icol_charge}];
    res{i,2} = temp_pep;
    res{i,3} = segment{column_idxs.icol_charge};
    res{i,4} = segment{column_idxs.icol_quant};
end

% Unique the name of the peptide list
unique_name = unique(res(:,1));

% Sum the quantification of the same peptide
uniqued_res = cell(length(unique_name),4);
for i = 1:length(unique_name)
    % Find the same peptide with the same charge
    idx_same_pep = strcmp(res(:,1),unique_name{i});

    % The name (sequence+charge), sequence, charge, quantification (sum)
    uniqued_res{i,1} = unique_name{i};
    uniqued_res{i,2} = res{find(idx_same_pep,1),2};
    uniqued_res{i,3} = res{find(idx_same_pep,1),3};
    uniqued_res{i,4} = sum(str2double({res{idx_same_pep,4}})); %#ok<CCAT1>
end

% The final result
res = uniqued_res;
end



function res = read_site_file(res_path, column_idxs)
% Read the site-level result file
% Input:
%   res_path
%       the path of the site-level result file
% Output:
%   res
%       the site-level result
%       [site_name, quantification, peptides(cell)]

% Open the file
fin = fopen(res_path,'r');
if fin == -1
    error(['Cannot open the result file of site summary: "', res_path, '"!']);
end

% Read the file and check how many sites are there
site_nums = 0;
strline = fgetl(fin);
while ~isequal(strline,-1)
    if ~isequal(strline(1),'*')
        site_nums = site_nums + 1;
    end
    strline = fgetl(fin);
end

% Rewind the file
frewind(fin);

% Initialize the table
res = cell(site_nums,3);

% Read the file and record the table
site_idx = 0;
strline = fgetl(fin);
while ~isequal(strline,-1)
    if ~isequal(strline(1),'*')
        % The site-started line, create a new site
        site_idx = site_idx + 1;
        segment = split(strline);
        res{site_idx,1} = segment{column_idxs.icol_site};   % The site name
        res{site_idx,2} = segment{column_idxs.icol_pep};    % The quantification
        res{site_idx,3} = {};
    else
        % The peptide line, record the peptide
        res{site_idx,3} = [res{site_idx,3}; strline];
    end
    strline = fgetl(fin);
end

% Close the file
fclose(fin);
end



function res = str_cmp_dict(str1, str2)
% compare two string using dictionary order
% Input:
%   str1, str2
%       two input strings
% Output:
%   res
%       0 when two input strings are the same
%       negative when the first is smaller in dictionary order
%       positive when the second is smaller in dictionary order
len1 = length(str1);
len2 = length(str2);
len_min = min([len1,len2]);
for i = 1:len_min
    if str1(i)-str2(i) ~= 0
        res = str1(i)-str2(i);
        return;
    end
end
res = len1-len2;
return;
end