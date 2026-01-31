function write_result(obj)
% Write the result to the file
% Input:
%   obj (CMergePairs)
%       merge instance with m_result populated

% Open the file
fout = fopen(obj.m_output_path,'w');
if fout == -1
    error(['Cannot open the output file: "', obj.m_output_path, '"!']);
end

% Write the header
fprintf(fout, 'Sites\tSequence\tCharge');
for idx_gt = 1:size(obj.m_group_titles, 1)
    fprintf(fout, '\t%s\t%s\tRatio\tNormalized ratio', obj.m_group_titles{idx_gt,1}, obj.m_group_titles{idx_gt,2});
end
fprintf(fout, '\n');

% Write the result to file
for idx_site = 1:size(obj.m_result,1)
    % Write the site name
    fprintf(fout, '%s', obj.m_result{idx_site,1});
    sum_pairs = zeros(size(obj.m_group_titles));
    
    % Write each peptide belong to the site
    pep_charge_quant_temp = obj.m_result{idx_site,2};
    for idx_pcqt = 1:size(pep_charge_quant_temp,1)
        % Write the peptide sequence and charge
        fprintf(fout, '\t%s\t%s', pep_charge_quant_temp{idx_pcqt,1}, pep_charge_quant_temp{idx_pcqt,2});
        for idx_pair = 1:size(pep_charge_quant_temp,2)-2    % Start from the third column
            quant_pair_temp = pep_charge_quant_temp{idx_pcqt, idx_pair+2};
            % Write the quantification pair (skip the empty pair), and add two tab after each pair 
            %   for manual calculating the ratio and normalized ratio
            if isempty(quant_pair_temp)
                fprintf(fout, '\t\t\t\t');
            else
                fprintf(fout, '\t%s\t%s\t\t', quant_pair_temp{1}, quant_pair_temp{2});
                sum_pairs(idx_pair, 1) = sum_pairs(idx_pair, 1) + str2double(quant_pair_temp{1});
                sum_pairs(idx_pair, 2) = sum_pairs(idx_pair, 2) + str2double(quant_pair_temp{2});
            end
        end
        fprintf(fout, '\n');
    end

    % Write the sum of each group
    fprintf(fout, '\tSUM\t');
    for idx_gt = 1:size(obj.m_group_titles, 1)
        % Comparison of ach pair should promise that both of them are empty or not in the same time
        if sum_pairs(idx_gt, 1) == 0
            fprintf(fout, '\t\t\t\t');
        else
            % Write the sum of each group (skip the empty pair) and the ratio, add one tab after each pair
            %   for manual calculating the normalized ratio
            fprintf(fout, '\t%f\t%f\t%f\t', sum_pairs(idx_gt, 1), sum_pairs(idx_gt, 2), ...
                sum_pairs(idx_gt, 1)/sum_pairs(idx_gt, 2));
        end
    end
    fprintf(fout, '\n');
end

% Close the file
fclose(fout);
end

