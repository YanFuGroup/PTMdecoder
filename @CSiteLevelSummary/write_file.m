function write_file(obj)
% Write to file
% Input:
%   obj (CSiteLevelSummary)
%       site-level summarizer instance

% Interested site-level
keys_interested = keys(obj.m_result_output_index);
fout = fopen(obj.m_output_path_interested, 'w');
if fout < 0
    error(['The file "', obj.m_output_path_interested, '" is not available.']);
end
for i_key = 1:length(keys_interested)
    fprintf(fout, '%s\t%f\n', keys_interested{i_key},...
        obj.m_result_output_sum{obj.m_result_output_index(keys_interested{i_key})});
    fprintf(fout, obj.m_result_output_string{obj.m_result_output_index(keys_interested{i_key})});
end
fclose(fout);

% Uninterested site-level
fout = fopen(obj.m_output_path_uninterested, 'w');
if fout < 0
    error(['The file "', obj.m_output_path_uninterested, '" is not available.']);
end
fclose(fout);
writecell(obj.m_result_uninterested_string, obj.m_output_path_uninterested,...
    'FileType', 'text', 'QuoteStrings', false);
end