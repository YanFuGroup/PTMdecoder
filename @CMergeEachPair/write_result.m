function write_result(obj)
% Write the merged result to file

% Open the file
fout = fopen(obj.m_output_path,'w');
if fout == -1
    error(['Cannot open the output file: "', obj.m_output_path, '"!']);
end

% Write the result to file
fprintf(fout, 'Sites\tSequence\tCharge\t%s\t%s\n', obj.m_group_titles{1}, obj.m_group_titles{2});
for i = 1:size(obj.m_result,1)
    % Write the site name
    fprintf(fout, '%s', obj.m_result{i,1});

    % Write the quantification of each peptide
    for j = 1:size(obj.m_result{i,2},1)
        fprintf(fout, '\t%s\t%s\t%f\t%f\n', obj.m_result{i,2}{j,1}, obj.m_result{i,2}{j,2},...
            obj.m_result{i,2}{j,3}, obj.m_result{i,2}{j,4});
    end
end

% Close the file
fclose(fout);
end